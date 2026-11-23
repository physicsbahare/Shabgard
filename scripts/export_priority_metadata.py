#!/usr/bin/env python3
"""Exact, archive-footprint-only priority metadata for Shabgard.

No substring-search of arbitrary CSV cells: program, mission and instrument
must be joined using the observation-level obsid. In particular, a file size
containing '7814' is NOT evidence for JWST program 7814.
"""
import argparse
import collections
import csv
import json
from pathlib import Path

GROUPS = ("JWST_7814_NIRCAM", "MIRI_IMAGE", "HST_WFC3_IR", "JWST_1727_NIRCAM")


def group_for_observation(row):
    if row.get("archive_footprint") != "inside_archive_footprint":
        return None
    mission = row.get("mission", "").upper()
    inst = row.get("instrument", "").upper()
    proposal = row.get("proposal", "").strip()
    kind = row.get("dataproduct_type", "").lower()
    if kind and kind != "image":
        return None
    if mission == "JWST" and "NIRCAM" in inst and "GRISM" not in inst:
        if proposal == "7814":
            return "JWST_7814_NIRCAM"
        if proposal == "1727":
            return "JWST_1727_NIRCAM"
    if mission == "JWST" and "MIRI/IMAGE" in inst:
        return "MIRI_IMAGE"
    if mission == "HST" and "WFC3/IR" in inst:
        return "HST_WFC3_IR"
    return None


def is_valid_product(row, group):
    fname = row.get("filename", "").lower()
    if row.get("rights", "PUBLIC").upper() not in ("PUBLIC", ""):
        return False
    if row.get("product_type", "").upper() not in ("SCIENCE", ""):
        return False
    if group in ("JWST_7814_NIRCAM", "JWST_1727_NIRCAM"):
        proposal = "07814" if group == "JWST_7814_NIRCAM" else "01727"
        return fname.startswith("jw" + proposal) and fname.endswith("_cal.fits") and "_nrc" in fname
    if group == "MIRI_IMAGE":
        return fname.endswith("_cal.fits") and "_mirimage_" in fname
    if group == "HST_WFC3_IR":
        # Only individual WFC3 exposures count as detector-level candidates.
        # HLA DRZ products can be retained for separate mosaic checks, NOT here.
        return fname.endswith(("_flt.fits", "_flc.fits"))
    return False


def load_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, columns):
    with path.open("w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        wr.writeheader()
        wr.writerows(rows)


def export(root, out):
    observations = load_csv(root / "observations.csv")
    products = load_csv(root / "science_products.csv")
    by_obsid = collections.defaultdict(list)
    selected_obs = []
    for row in observations:
        group = group_for_observation(row)
        if not group:
            continue
        obsid = row.get("obsid", "")
        if not obsid:
            continue
        selected_obs.append({"priority_group": group, **row})
        by_obsid[obsid].append((group, row))
    selected_prod = []
    counts = collections.Counter()
    unmatched = collections.Counter()
    for product in products:
        oid = product.get("query_obsid", "")
        if not oid or oid not in by_obsid:
            if product.get("filename", "").lower().endswith(("_cal.fits", "_flt.fits", "_flc.fits")):
                unmatched["calibrated_product_without_selected_obsid_link"] += 1
            continue
        for group, obs in by_obsid[oid]:
            if not is_valid_product(product, group):
                continue
            selected_prod.append({
                "priority_group": group,
                "proposal": obs.get("proposal", ""),
                "instrument": obs.get("instrument", ""),
                "filter": obs.get("filter", ""),
                "date_start_utc": obs.get("date_start_utc", ""),
                "observation_obsid": oid,
                **product,
            })
            counts[group] += 1
    # Product names alone don't prove detector coverage; avoid duplicate rows
    # if one obsid occurs multiple times with identical group/filter metadata.
    selected_prod = list({(r["priority_group"], r["data_uri"], r["filter"]): r
                          for r in selected_prod}.values())
    selected_obs.sort(key=lambda r: (GROUPS.index(r["priority_group"]),
                                     r.get("date_start_utc", ""), r.get("obsid", "")))
    selected_prod.sort(key=lambda r: (GROUPS.index(r["priority_group"]),
                                      r.get("date_start_utc", ""), r.get("filter", ""),
                                      r.get("filename", "")))
    out.mkdir(parents=True, exist_ok=True)
    obs_columns = ["priority_group"] + list(observations[0]) if observations else ["priority_group"]
    prod_columns = ["priority_group", "proposal", "instrument", "filter",
                    "date_start_utc", "observation_obsid"] + (list(products[0]) if products else [])
    write_csv(out / "priority_observations.csv", selected_obs, obs_columns)
    write_csv(out / "priority_calibrated_products.csv", selected_prod, prod_columns)
    for group in GROUPS:
        write_csv(out / (group.lower() + "_observations.csv"),
                  [r for r in selected_obs if r["priority_group"] == group], obs_columns)
        write_csv(out / (group.lower() + "_products.csv"),
                  [r for r in selected_prod if r["priority_group"] == group], prod_columns)
    summary = {
        "method": "Exact observation metadata group selection, then query_obsid join",
        "observations_total": len(observations),
        "science_products_total": len(products),
        "priority_observations": dict(collections.Counter(r["priority_group"] for r in selected_obs)),
        "calibrated_products": dict(collections.Counter(r["priority_group"] for r in selected_prod)),
        "unmatched": dict(unmatched),
        "scientific_status": "ARCHIVE FOOTPRINT METADATA ONLY; no calibrated detector coverage or detection",
        "warning": "Some products lack query_obsid; they are deliberately not assumed to overlap Shabgard."
    }
    (out / "schema_and_counts.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--results", default="results")
    p.add_argument("--out", default="priority_metadata")
    a = p.parse_args()
    print(json.dumps(export(Path(a.results), Path(a.out)), indent=2))


if __name__ == "__main__":
    main()
