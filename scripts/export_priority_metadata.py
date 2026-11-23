#!/usr/bin/env python3
"""Generate small, inspectable priority subsets from the large MAST CSV outputs.

This is metadata triage only. It does not establish detector coverage or detection.
"""
import argparse
import csv
import json
from pathlib import Path

PRIORITIES = [
    ("JWST_7814", ("7814",)),
    ("MIRI", ("MIRI", "F770W", "F1000W", "F1280W", "F1500W", "F1800W", "F2100W", "F2550W")),
    ("HST_WFC3_IR", ("WFC3/IR", "WFC3_IR", "WFC3 IR", "F125W", "F140W", "F160W")),
]
MAX_ROWS = 80


def read_and_select(path, output_dir):
    if not path.exists():
        return {"file": str(path), "error": "missing"}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return {"file": str(path), "error": "no CSV header"}
        columns = reader.fieldnames
        selected = {name: [] for name, _ in PRIORITIES}
        totals = {name: 0 for name, _ in PRIORITIES}
        for row in reader:
            # Inspect all metadata fields, not just an assumed column schema.
            value = " ".join(str(v or "") for v in row.values()).upper()
            for name, needles in PRIORITIES:
                if any(needle.upper() in value for needle in needles):
                    totals[name] += 1
                    if len(selected[name]) < MAX_ROWS:
                        selected[name].append(row)
    for name, rows in selected.items():
        target = output_dir / f"{path.stem}_{name.lower()}_sample.csv"
        with target.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
    return {"file": str(path), "columns": columns,
            "matching_rows": totals, "sample_limit_per_priority": MAX_ROWS}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", default="results")
    ap.add_argument("--out", default="results/priority_metadata")
    args = ap.parse_args()
    directory = Path(args.out)
    directory.mkdir(parents=True, exist_ok=True)
    reports = [read_and_select(Path(args.results) / filename, directory)
               for filename in ("observations.csv", "science_products.csv")]
    (directory / "schema_and_counts.json").write_text(
        json.dumps(reports, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    for report in reports:
        print(report["file"], report.get("matching_rows", report.get("error")))


if __name__ == "__main__":
    main()
