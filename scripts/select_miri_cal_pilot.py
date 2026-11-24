#!/usr/bin/env python3
"""Select a bounded, auditable MIRI CAL pilot from MAST, not mosaics.

Priority F770W, F1000W, F1280W, F1500W. F770W metadata is queried
directly when absent from the persisted product table. Every selected
CAL is a geometric *candidate* until its individual gWCS/SCI/DQ is read.
"""
import argparse
import csv
import json
from pathlib import Path

OBS = [
    ("F770W", "109374093", "1837", "2022-12-26"),
    ("F770W", "200464999", "1727", "2024-01-04"),
    ("F770W", "200763671", "1727", "2024-01-06"),
    ("F1000W", "234492195", "5893", "2024-12-21"),
    ("F1280W", "372523982", "7814", "2025-12-29"),
    ("F1500W", "372524111", "7814", "2025-12-29"),
]


def value(row, field):
    try:
        v = str(row[field])
        return "" if v in ("None", "nan", "--") else v
    except (KeyError, IndexError):
        return ""


def archive_products(obsid):
    from astroquery.mast import Observations
    rows = []
    for p in Observations.get_product_list(obsid):
        name = value(p, "productFilename")
        if name.lower().endswith("_mirimage_cal.fits"):
            rows.append(dict(filename=name, data_uri=value(p, "dataURI"),
                             rights=value(p, "dataRights"), size_bytes=value(p, "size"),
                             subgroup=value(p, "productSubGroupDescription")))
    return rows


def select(root, out, max_per_filter, include_query, filters=None, pointings=None, obsids=None):
    existing = {}
    path = root / "science_products.csv"
    if path.exists():
        with path.open(newline="") as f:
            for p in csv.DictReader(f):
                if p.get("filename", "").lower().endswith("_mirimage_cal.fits"):
                    existing.setdefault(p.get("query_obsid", ""), []).append(p)
    rows, report = [], []
    for filt, obsid, proposal, date in OBS:
        if filters is not None and filt not in filters:
            continue
        if obsids is not None and obsid not in obsids:
            continue
        products = existing.get(obsid, [])
        source = "persisted_science_products"
        if not products and include_query:
            try:
                products = archive_products(obsid)
                source = "live_MAST_product_query"
            except Exception as exc:
                source = "live_MAST_failed:" + str(exc)[:150]
        products = [p for p in products if p.get("filename", "").endswith("_mirimage_cal.fits")
                    and p.get("rights", "PUBLIC").upper() in ("PUBLIC", "")
                    and (pointings is None or p["filename"].split("_")[0] in pointings)]
        # One first exposure from every distinct pointing, then remaining
        # dithers. This covers distinct sky tiles earlier than serial order.
        by_pointing = {}
        for p in sorted(products, key=lambda p: p["filename"]):
            prefix = p["filename"].split("_")[0]
            by_pointing.setdefault(prefix, []).append(p)
        candidates = [xs[0] for xs in by_pointing.values()]
        candidates += [p for xs in by_pointing.values() for p in xs[1:]]
        # Sample evenly across pointings; do not treat sampled noncoverage
        # as proof that the entire observation missed the target.
        n = min(max_per_filter, len(candidates))
        selected = []
        if n:
            start = min(n, len(by_pointing))
            selected = [candidates[int(i * len(by_pointing) / start)] for i in range(start)]
            if start < n:
                selected.extend(candidates[len(by_pointing):len(by_pointing) + n - start])
        for p in selected:
            rows.append(dict(filter=filt, obsid=obsid, proposal=proposal,
                             date_utc=date, filename=p["filename"],
                             data_uri=p.get("data_uri") or "mast:JWST/product/" + p["filename"],
                             size_bytes=p.get("size_bytes", ""), source=source))
        report.append(dict(filter=filt, obsid=obsid, date_utc=date, source=source,
                           total_cal=len(products), selected=len(selected),
                           sampled_fraction=(len(selected) / len(products) if products else 0)))
    out.mkdir(parents=True, exist_ok=True)
    with (out / "miri_cal_selection.csv").open("w", newline="") as f:
        columns = ["filter", "obsid", "proposal", "date_utc", "filename",
                   "data_uri", "size_bytes", "source"]
        w = csv.DictWriter(f, fieldnames=columns)
        w.writeheader()
        w.writerows(rows)
    (out / "miri_cal_selection_summary.json").write_text(
        json.dumps({"observations": report,
                    "caveat": "Sampling CAL candidates; negative coverage in sampled files cannot rule out target in untested pointings"},
                   indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if not rows:
        raise RuntimeError("No MIRI CAL products found; investigate archive query")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--out", default="miri_results")
    ap.add_argument("--max-per-observation", type=int, default=4)
    ap.add_argument("--no-live-query", action="store_true")
    ap.add_argument("--filters", nargs="+", choices=("F770W","F1000W","F1280W","F1500W"), default=None)
    ap.add_argument("--pointings", nargs="+", default=None, help="Exact JWST visit+pointing prefixes from independently tested CAL")
    ap.add_argument("--obsids", nargs="+", default=None, help="Exact MAST observation IDs")
    args = ap.parse_args()
    select(Path(args.results), Path(args.out), args.max_per_observation,
           not args.no_live_query, args.filters, args.pointings, args.obsids)


if __name__ == "__main__":
    main()
