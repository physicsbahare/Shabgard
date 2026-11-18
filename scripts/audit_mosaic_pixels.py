#!/usr/bin/env python3
"""Read small FITS SCI/WHT/DQ cutouts remotely; distinguish WCS from usable pixels.
A positive result is mosaic pixel coverage only, NOT a detector exposure or detection.
"""
import argparse
import csv
import json
import math
import signal
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote
import numpy as np
from astropy.io import fits

FIELDS=["data_uri","product","mission","instrument","filter","date_utc",
        "x_motion_proxy","y_motion_proxy","pixel_status","sci_finite_pixels",
        "sci_total_pixels","sci_median","sci_mad","wht_positive_pixels",
        "dq_good_pixels","pixel_cutout_size","error"]
def alarm(sig, frame):
    raise TimeoutError("Remote pixel range-read timed out")
def sample(row, timeout):
    uri=row["data_uri"]
    url="https://mast.stsci.edu/api/v0.1/Download/file?uri="+quote(uri,safe=":/")
    x=float(row["x_motion_proxy"]);y=float(row["y_motion_proxy"])
    signal.alarm(timeout)
    try:
        with fits.open(url,use_fsspec=True,
                       fsspec_kwargs={"block_size":262144},
                       lazy_load_hdus=True,memmap=False) as hdul:
            sci=next((h for h in hdul if h.name=="SCI" and h.header.get("NAXIS",0)>=2),None)
            # HLA drizzled WFC3/IR products may store their science image
            # in the primary HDU, unlike JWST i2d products.
            if sci is None and hdul[0].header.get("NAXIS",0)>=2:
                sci=hdul[0]
            if sci is None: return dict(pixel_status="no_image_hdu",error="No 2D SCI or primary image")
            nx=int(sci.header["NAXIS1"]);ny=int(sci.header["NAXIS2"])
            ix=int(round(x));iy=int(round(y))
            if ix<0 or iy<0 or ix>=nx or iy>=ny:
                return dict(pixel_status="outside_sci_array",error="")
            radius=5
            x0=max(0,ix-radius);x1=min(nx,ix+radius+1)
            y0=max(0,iy-radius);y1=min(ny,iy+radius+1)
            # FITS section requests only the local pixels over HTTP ranges.
            a=np.asarray(sci.section[y0:y1,x0:x1],dtype=float)
            valid=np.isfinite(a)
            n=int(valid.sum());total=int(a.size)
            out=dict(sci_finite_pixels=n,sci_total_pixels=total,
                     pixel_cutout_size=f"{a.shape[1]}x{a.shape[0]}",
                     sci_median=float(np.median(a[valid])) if n else "",
                     sci_mad=float(np.median(abs(a[valid]-np.median(a[valid])))) if n else "")
            for ext,field,rule in (("WHT","wht_positive_pixels",lambda z: np.isfinite(z)&(z>0)),
                                   ("DQ","dq_good_pixels",lambda z: np.isfinite(z)&(z==0))):
                h=next((h for h in hdul if h.name==ext and h.header.get("NAXIS",0)>=2),None)
                if h is not None and h.header.get("NAXIS1")==nx and h.header.get("NAXIS2")==ny:
                    b=np.asarray(h.section[y0:y1,x0:x1])
                    out[field]=int(rule(b).sum())
            if n==0: out["pixel_status"]="no_finite_sci_pixels"
            elif out.get("wht_positive_pixels",-1)==0: out["pixel_status"]="zero_weight"
            else: out["pixel_status"]="finite_mosaic_sci_pixels"
            out["error"]=""
            return out
    finally: signal.alarm(0)
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",default="results")
    p.add_argument("--limit",type=int,default=12)
    p.add_argument("--timeout",type=int,default=30)
    args=p.parse_args()
    root=Path(args.out)
    candidates=list(csv.DictReader((root/"mosaic_wcs_audit.csv").open()))
    path=root/"mosaic_pixel_audit.csv"
    old=list(csv.DictReader(path.open())) if path.exists() else []
    existing={r["data_uri"]:r for r in old}
    signal.signal(signal.SIGALRM,alarm)
    n=0
    for r in candidates:
        if r.get("wcs_status")!="inside_mosaic_image_wcs":continue
        if r["data_uri"] in existing and existing[r["data_uri"]].get("pixel_status") not in ("no_sci_hdu","no_image_hdu","remote_error_or_timeout"):continue
        out={k:r.get(k,"") for k in FIELDS}
        try:out.update(sample(r,args.timeout))
        except Exception as exc:out.update(pixel_status="remote_error_or_timeout",error=(type(exc).__name__+": "+str(exc))[:250])
        existing[r["data_uri"]]=out
        n+=1
        if n>=args.limit:break
    rows=list(existing.values())
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,extrasaction="ignore")
        w.writeheader();w.writerows(rows)
    counts={}
    for r in rows:counts[r["pixel_status"]]=counts.get(r["pixel_status"],0)+1
    report={"updated_utc":datetime.now(timezone.utc).isoformat(),
            "total_products_sampled":len(rows),"sampled_this_run":n,
            "pixel_status_counts":counts,
            "note":"Finite SCI mosaic pixels do not prove detector exposure coverage or a target detection; WHT/DQ checks where available. No SNR/centroid claim."}
    (root/"mosaic_pixel_summary.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
if __name__=="__main__":main()
