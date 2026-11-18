#!/usr/bin/env python3
"""Incremental MAST product metadata fetch, with persistent cursor and retry.

Reads already-vetted archive-footprint records; does not assert detector coverage.
A small number of obsids is queried per run to avoid MAST timeouts.
"""
import argparse
import csv
import json
import time
import signal
from datetime import datetime, timezone
from pathlib import Path
from astroquery.mast import Observations

FIELDS = ["obsid", "filename", "data_uri", "product_type", "subgroup", "size_bytes", "rights"]
SUFFIX = ("_cal.fits", "_rate.fits", "_i2d.fits", "_drz.fits", "_drc.fits",
          "_flt.fits", "_flc.fits", "_s2d.fits", "_x1d.fits", "_c0f.fits")
def value(row, name):
    try:
        x = row[name]
        return "" if x is None or str(x) in ("--", "nan", "None") else str(x)
    except (KeyError, IndexError, ValueError):
        return ""

def write_csv(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results")
    ap.add_argument("--batch-size", type=int, default=3)
    args = ap.parse_args()
    root = Path(args.out)
    observations = list(csv.DictReader((root / "observations.csv").open()))
    obsids = sorted({r["obsid"] for r in observations
                     if r.get("archive_footprint") == "inside_archive_footprint"
                     and r.get("obsid", "").isdigit()})
    cursor_file = root / "product_progress.json"
    progress = json.loads(cursor_file.read_text()) if cursor_file.exists() else {}
    completed = set(progress.get("completed_obsids", []))
    failures = dict(progress.get("failed_obsids", {}))
    # Fast-science priority: new JWST/NIRCam 2026, all MIRI, other
    # JWST/NIRCam, then HST/WFC3-IR, then all other archival records.
    # These are metadata candidates, not verified source detections.
    by_id = {r["obsid"]: r for r in observations
             if r.get("archive_footprint") == "inside_archive_footprint"}
    def rank(oid):
        r = by_id.get(oid, {})
        mission, instrument = r.get("mission", ""), r.get("instrument", "").upper()
        proposal = r.get("proposal", "")
        if mission == "JWST" and proposal == "7814" and "NIRCAM" in instrument: return 0
        if mission == "JWST" and "MIRI" in instrument: return 1
        if mission == "JWST" and "NIRCAM" in instrument: return 2
        if mission == "HST" and "WFC3/IR" in instrument: return 3
        return 4
    pending = sorted([o for o in obsids if o not in completed],
                     key=lambda oid: (rank(oid), oid))
    def expire(signum, frame):
        raise TimeoutError("MAST product lookup exceeded 60 seconds for an obsid")
    signal.signal(signal.SIGALRM, expire)
    product_path = root / "science_products.csv"
    old = list(csv.DictReader(product_path.open())) if product_path.exists() else []
    indexed = {r["data_uri"] or r["filename"]: r for r in old}
    done_this_run = []
    for oid in pending[:max(1, args.batch_size)]:
        try:
            # One obsid at a time, so failures cannot discard other batches.
            signal.alarm(60)
            pp = Observations.get_product_list(oid)
            for p in pp:
                filename = value(p, "productFilename")
                if not filename.lower().endswith(SUFFIX):
                    continue
                record = dict(obsid=value(p, "obsID"), filename=filename,
                              data_uri=value(p, "dataURI"),
                              product_type=value(p, "productType"),
                              subgroup=value(p, "productSubGroupDescription"),
                              size_bytes=value(p, "size"),
                              rights=value(p, "dataRights"))
                indexed[record["data_uri"] or filename] = record
            completed.add(oid)
            failures.pop(oid, None)
            done_this_run.append(oid)
        except Exception as e:
            failures[oid] = str(e)[:300]
        finally:
            signal.alarm(0)
        time.sleep(0.25)
    write_csv(product_path, sorted(indexed.values(), key=lambda r:(r["obsid"],r["filename"])))
    progress = dict(updated_utc=datetime.now(timezone.utc).isoformat(),
                    source="MAST product metadata", total_candidate_obsids=len(obsids),
                    completed_obsids=sorted(completed), failed_obsids=failures,
                    remaining=len([x for x in obsids if x not in completed]),
                    science_product_files=len(indexed),
                    detector_coverage_verified=False,
                    source_detected=False,
                    processed_this_run=done_this_run,
                    priority="JWST 7814 NIRCam, JWST MIRI, other JWST NIRCam, HST WFC3/IR, other")
    cursor_file.write_text(json.dumps(progress, indent=2)+"\n")
    print(json.dumps({k:v for k,v in progress.items() if k != "completed_obsids"}, indent=2))

if __name__ == "__main__":
    main()
