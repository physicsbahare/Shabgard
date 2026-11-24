#!/usr/bin/env python3
"""Download only selected public MAST MIRI CAL files with a checked manifest.

Avoids downloading full archive or mistaking a failed HTTP response for FITS.
Each file is <=~30 MB in currently indexed MIRI program selections.
"""
import argparse
import csv
import json
import time
from pathlib import Path
from urllib.parse import quote

import requests


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--selection", default="miri_results/miri_cal_selection.csv")
    p.add_argument("--out", default="miri_data")
    p.add_argument("--report", default="miri_results/miri_download_manifest.csv")
    p.add_argument("--max-files", type=int, default=24)
    args = p.parse_args()
    root = Path(args.out)
    root.mkdir(parents=True, exist_ok=True)
    with open(args.selection, newline="") as f:
        targets = list(csv.DictReader(f))[:args.max_files]
    if not targets:
        raise SystemExit("No selected MIRI products to download")
    reports = []
    session = requests.Session()
    for r in targets:
        name = r["filename"]
        uri = r["data_uri"]
        outfile = root / name
        rec = {"filter": r["filter"], "filename": name, "data_uri": uri,
               "status": "error", "bytes": 0, "error": ""}
        if "/" in name or not name.endswith("_mirimage_cal.fits"):
            rec["error"] = "Disallowed filename"
            reports.append(rec)
            continue
        url = "https://mast.stsci.edu/api/v0.1/Download/file?uri=" + quote(uri, safe=":/")
        for attempt in range(3):
            try:
                if outfile.exists() and outfile.stat().st_size > 1_000_000:
                    with outfile.open("rb") as f:
                        isfits = f.read(24).startswith(b"SIMPLE  =")
                    if isfits:
                        rec.update(status="ok_cached", bytes=outfile.stat().st_size)
                        break
                with session.get(url, stream=True, timeout=(20, 170)) as rsp:
                    rsp.raise_for_status()
                    nbytes = 0
                    part = outfile.with_suffix(".download")
                    with part.open("wb") as out:
                        for data in rsp.iter_content(chunk_size=1 << 20):
                            if data:
                                out.write(data)
                                nbytes += len(data)
                    if nbytes < 1_000_000 or not part.open("rb").read(24).startswith(b"SIMPLE  ="):
                        part.unlink(missing_ok=True)
                        raise ValueError("Response is not a complete FITS CAL file")
                    part.replace(outfile)
                    rec.update(status="ok", bytes=nbytes)
                    break
            except Exception as exc:
                rec["error"] = type(exc).__name__ + ": " + str(exc)[:180]
                if attempt < 2:
                    time.sleep(3 * (attempt + 1))
        reports.append(rec)
        print(rec, flush=True)
    dest = Path(args.report)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(reports[0]))
        w.writeheader()
        w.writerows(reports)
    print(json.dumps({"requested": len(reports),
                      "downloaded": sum(r["status"].startswith("ok") for r in reports),
                      "failed": sum(not r["status"].startswith("ok") for r in reports)}, indent=2))
    if not any(r["status"].startswith("ok") for r in reports):
        raise SystemExit("No MIRI CAL files downloaded; inspect MAST errors")


if __name__ == "__main__":
    main()
