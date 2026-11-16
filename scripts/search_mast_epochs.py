#!/usr/bin/env python3
"""Inventory HST and JWST observation/SCIENCE product metadata near Shabgard.

The intersection with archive s_region is an OBSERVATION-FOOTPRINT test, not
a guarantee that the source lies in a particular calibrated detector's pixels.
"""
import argparse
import csv
import json
import math
import re
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from astroquery.mast import Observations
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
from shapely.geometry import Point, Polygon

SCIENCE_SUFFIXES = (
    "_cal.fits", "_rate.fits", "_i2d.fits", "_drz.fits", "_drc.fits",
    "_flt.fits", "_flc.fits", "_s2d.fits", "_x1d.fits", "_c0f.fits",
)
MAIN = ("JWST", "HST", "HLA")


def val(x):
    if x is None or str(x) in ("--", "nan", "None"):
        return ""
    return str(x)


def field(row, name):
    try:
        return val(row[name])
    except (KeyError, ValueError, IndexError):
        return ""


def fnum(x):
    try:
        n = float(x)
        return n if math.isfinite(n) else None
    except (TypeError, ValueError):
        return None


def xy(ra, dec, target_ra, target_dec):
    return (((ra - target_ra + 180) % 360 - 180) * math.cos(math.radians(target_dec)),
            dec - target_dec)


def footprint_status(region, ra, dec):
    """Inside/outside for simple STC-S polygons/circles; unknown otherwise."""
    if not region:
        return "unknown_no_s_region"
    r = region.upper()
    point = Point(0, 0)
    found_shape = False
    try:
        for seg in re.split(r"\bPOLYGON\b", r)[1:]:
            seg = re.split(r"\b(?:CIRCLE|BOX|UNION|INTERSECTION|DIFFERENCE)\b|\)", seg)[0]
            nums = [float(z) for z in re.findall(r"[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[Ee][-+]?\d+)?", seg)]
            if len(nums) < 6 or len(nums) % 2:
                continue
            pts = [xy(nums[i], nums[i + 1], ra, dec) for i in range(0, len(nums), 2)]
            poly = Polygon(pts)
            if not poly.is_valid:
                poly = poly.buffer(0)
            if not poly.is_empty:
                found_shape = True
                if poly.covers(point) or poly.distance(point) < 1e-6:
                    return "inside_archive_footprint"
        for seg in re.split(r"\bCIRCLE\b", r)[1:]:
            seg = re.split(r"\b(?:POLYGON|BOX|UNION|INTERSECTION|DIFFERENCE)\b|\)", seg)[0]
            nums = [float(z) for z in re.findall(r"[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[Ee][-+]?\d+)?", seg)]
            if len(nums) >= 3:
                found_shape = True
                xx, yy = xy(nums[0], nums[1], ra, dec)
                if math.hypot(xx, yy) <= abs(nums[2]) + 1e-6:
                    return "inside_archive_footprint"
        return "outside_archive_footprint" if found_shape else "unknown_region_syntax"
    except Exception:
        return "unknown_region_parse_error"


def mjd_iso(value):
    v = fnum(value)
    if v is None:
        return ""
    try:
        return Time(v, format="mjd", scale="utc").isot
    except Exception:
        return ""


def write_csv(path, records, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        w.writerows(records)


def query_with_retry(position, radius_arcmin):
    last = None
    for attempt in range(3):
        try:
            return Observations.query_region(position, radius=radius_arcmin * u.arcmin)
        except Exception as e:
            last = e
            print(f"MAST query attempt {attempt+1} failed: {e}", file=sys.stderr)
            time.sleep(4 * (attempt + 1))
    raise RuntimeError(f"MAST observation query failed: {last}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="target.json")
    ap.add_argument("--out", default="results")
    args = ap.parse_args()
    conf = json.loads(Path(args.target).read_text())
    ra, dec = float(conf["ra_deg"]), float(conf["dec_deg"])
    radius = float(conf.get("search_radius_arcmin", 4))
    root = Path(args.out)
    root.mkdir(parents=True, exist_ok=True)
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    status = dict(query_started_utc=now, ra_deg=ra, dec_deg=dec,
                  coordinate_status=conf.get("coordinate_status", "unknown"),
                  radius_arcmin=radius, source="MAST / astroquery.mast",
                  ok=False, product_query_complete=False)
    try:
        position = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
        all_rows = query_with_retry(position, radius)
        observations, others, eligible = [], [], []
        for i, row in enumerate(all_rows):
            mission = field(row, "obs_collection").upper()
            r = dict(
                mission=mission, instrument=field(row, "instrument_name"),
                filter=field(row, "filters"), proposal=field(row, "proposal_id"),
                obs_id=field(row, "obs_id"), obsid=field(row, "obsid"),
                date_start_utc=mjd_iso(field(row, "t_min")),
                date_end_utc=mjd_iso(field(row, "t_max")),
                t_min_mjd=field(row, "t_min"), t_max_mjd=field(row, "t_max"),
                exposure_s=field(row, "t_exptime"),
                obs_title=field(row, "obs_title"),
                calib_level=field(row, "calib_level"),
                dataproduct_type=field(row, "dataproduct_type"),
                rights=field(row, "dataRights"), s_region=field(row, "s_region"),
                archive_footprint=footprint_status(field(row, "s_region"), ra, dec),
            )
            if mission in MAIN:
                observations.append(r)
                if r["archive_footprint"] != "outside_archive_footprint":
                    eligible.append((i, r))
            else:
                others.append(r)
        obs_fields = list((observations or others or [dict(mission="")])[0].keys())
        write_csv(root / "observations.csv", observations, obs_fields)
        write_csv(root / "other_missions.csv", others, obs_fields)

        # Groups reflect archive observation start time, not exact individual
        # exposure MJD. Require exposure-level headers before parallax fitting.
        groups = defaultdict(list)
        for _, r in eligible:
            date = r["date_start_utc"][:10] or "unknown_date"
            key = (date, r["mission"], r["instrument"], r["filter"], r["proposal"])
            groups[key].append(r)
        epoch_records = []
        for (date, mission, inst, filt, proposal), records in sorted(groups.items()):
            epoch_records.append(dict(
                date_utc=date, mission=mission, instrument=inst, filter=filt,
                proposal=proposal, observation_records=len(records),
                archive_footprint_hits=sum(
                    r["archive_footprint"] == "inside_archive_footprint" for r in records),
                unknown_footprint=sum(
                    r["archive_footprint"].startswith("unknown") for r in records),
                obs_ids=";".join(sorted(set(r["obs_id"] for r in records))),
                assessment="CANDIDATE only: no detector/pixel or SNR validation",
            ))
        ep_fields = ["date_utc", "mission", "instrument", "filter", "proposal",
                     "observation_records", "archive_footprint_hits", "unknown_footprint",
                     "obs_ids", "assessment"]
        write_csv(root / "epoch_groups.csv", epoch_records, ep_fields)

        products, errors = [], []
        # Product metadata is required before per-exposure WCS and centroid work.
        # Request small batches so one timeout does not discard the whole inventory.
        for off in range(0, len(eligible), 30):
            indices = [i for i, _ in eligible[off:off + 30]]
            selected = all_rows[indices]
            try:
                pp = Observations.get_product_list(selected)
                for p in pp:
                    filename = field(p, "productFilename")
                    if not filename.lower().endswith(SCIENCE_SUFFIXES):
                        continue
                    products.append(dict(
                        obsid=field(p, "obsID"), filename=filename,
                        data_uri=field(p, "dataURI"),
                        product_type=field(p, "productType"),
                        subgroup=field(p, "productSubGroupDescription"),
                        size_bytes=field(p, "size"),
                        rights=field(p, "dataRights"),
                    ))
            except Exception as exc:
                errors.append(f"product batch starting {off}: {type(exc).__name__}: {exc}")
        prod_fields = ["obsid", "filename", "data_uri", "product_type",
                       "subgroup", "size_bytes", "rights"]
        # Product listing can return the same file from overlapping observation rows.
        uniq = {r["data_uri"] or r["filename"]: r for r in products}
        products = sorted(uniq.values(), key=lambda r: (r["obsid"], r["filename"]))
        write_csv(root / "science_products.csv", products, prod_fields)

        counts = Counter(r["mission"] for _, r in eligible)
        instrument_counts = Counter((r["mission"], r["instrument"]) for _, r in eligible)
        status.update(
            ok=True, query_completed_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
            raw_observation_count=len(all_rows),
            scoped_observation_count=len(observations),
            candidate_footprint_records=len(eligible),
            archive_footprint_hits=sum(r["archive_footprint"] == "inside_archive_footprint" for _, r in eligible),
            candidate_epoch_groups=len(epoch_records),
            science_product_files=len(products),
            product_query_complete=not errors,
            product_query_errors=errors,
        )
        lines = [
            "# Shabgard: MAST observation census",
            "", f"MAST query: {now} (UTC), radius {radius} arcmin",
            f"Coordinates: RA {ra:.8f} deg, Dec {dec:.8f} deg "
            f"(**{conf.get('coordinate_status', 'unverified')}**)",
            "",
            "> **Important:** These are archive-level footprint candidates, NOT verified",
            "> target detections or independent high-precision astrometric epochs.",
            "> The single-image detector WCS, gaps, SNR and target centroid remain to be checked.",
            "",
            f"- Total raw nearby archive observation records: **{len(all_rows)}**",
            f"- JWST/HST/HLA observation records: **{len(observations)}**",
            f"- Candidate overlap or unknown footprint: **{len(eligible)}**",
            f"- Explicit archive footprint hits: **{status['archive_footprint_hits']}**",
            f"- Date × mission × instrument × filter × program groups: **{len(epoch_records)}**",
            f"- Science FITS products indexed: **{len(products)}**",
            "",
            "## Mission/instrument breakdown",
            "",
        ]
        for (mission, inst), n in sorted(instrument_counts.items()):
            lines.append(f"- {mission} / {inst}: {n} candidate observation records")
        lines.extend(["", "## Candidate observing dates", "",
                      "| UTC date | Mission | Instrument | Filter | Program | Records |",
                      "|---|---|---|---|---|---:|"])
        for r in epoch_records:
            lines.append(f"| {r['date_utc']} | {r['mission']} | "
                         f"{r['instrument']} | {r['filter']} | {r['proposal']} | "
                         f"{r['observation_records']} |")
        lines.extend(["", "## Caveats and next checks", "",
                      "1. Confirm the target RA/Dec against the COSMOS2025 source row.",
                      "2. Verify each calibrated detector's exact pixel/gWCS footprint; MAST metadata polygons may span gaps.",
                      "3. Verify Shabgard is actually detected and centroidable in each relevant HST/MIRI/JWST filter.",
                      "4. Read per-exposure MJD and collapse dithers/same visits into independent astrophysical epochs.",
                      "5. Test parallax-factor spread, cross-filter systematics and Gaia/galaxy reference ties.",
                      ])
        if errors:
            lines.extend(["", "## Product-query partial failures", "",
                          "The observation census succeeded, but the product list is incomplete."])
            lines.extend(f"- {err}" for err in errors)
        (root / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    except Exception as e:
        status["error"] = f"{type(e).__name__}: {e}"
        (root / "REPORT.md").write_text(
            "# Shabgard archive query failed\n\n" + status["error"] +
            "\n\nNo new coverage can be inferred from this run.\n", encoding="utf-8")
        raise
    finally:
        (root / "status.json").write_text(
            json.dumps(status, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
