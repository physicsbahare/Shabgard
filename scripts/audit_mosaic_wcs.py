#!/usr/bin/env python3
"""FAST mosaic image-WCS screening; no claim of target detection or parallax."""
import argparse
import csv
import json
import math
import signal
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

FIELDS = ["mission","instrument","proposal","filter","obsid","obs_id","date_utc",
          "product","data_uri","wcs_status","wcs_hdu","x_catalog","y_catalog",
          "x_motion_proxy","y_motion_proxy","naxis1","naxis2",
          "pixel_scale_arcsec","catalog_inbounds","motion_proxy_inbounds",
          "quality_note","error"]
SUFFIXES = ("_i2d.fits","_drz.fits","_drc.fits")

def priority(obs):
    mission = obs.get("mission","")
    inst = obs.get("instrument","").upper()
    prop = obs.get("proposal","")
    if mission == "JWST" and prop == "7814" and "NIRCAM" in inst: return 0
    if mission == "JWST" and "MIRI" in inst: return 1
    if mission == "JWST" and "NIRCAM" in inst: return 2
    if mission == "HST" and "WFC3/IR" in inst: return 3
    return 4

def timeout_handler(sig, frame):
    raise TimeoutError("Remote FITS header request exceeded timeout")

def approximate_position(ra, dec, date):
    # Proper-motion F115W components Salmasi+2026 §5.3.2, NOT catalog epoch.
    # 2023.99 is only an explicit proxy; HST positions especially uncertain.
    try:
        delta = Time(date, format="isot", scale="utc").decimalyear - 2023.99
    except Exception:
        delta = 0.0
    return (ra - 70.36*delta/(3.6e6*math.cos(math.radians(dec))),
            dec - 53.89*delta/3.6e6)

def check(uri, ra, dec, obs, seconds):
    if not uri.startswith("mast:"):
        return dict(wcs_status="unknown_uri", error="No MAST URI")
    url = "https://mast.stsci.edu/api/v0.1/Download/file?uri=" + quote(uri,safe=":/")
    signal.alarm(seconds)
    try:
        # fsspec range requests only. Never call .data on full FITS image.
        with fits.open(url, use_fsspec=True,
                       fsspec_kwargs={"block_size":262144},
                       lazy_load_hdus=True, memmap=False) as hdus:
            for hdu in hdus:
                header = hdu.header
                if header.get("NAXIS",0) < 2 or "CTYPE1" not in header:
                    continue
                nx,ny = int(header.get("NAXIS1",0)),int(header.get("NAXIS2",0))
                if min(nx,ny) < 1: continue
                try:
                    w = WCS(header).celestial
                    x,y = w.all_world2pix(ra,dec,0)
                    mra,mdec = approximate_position(ra,dec,obs.get("date_start_utc",""))
                    xm,ym = w.all_world2pix(mra,mdec,0)
                    if not all(math.isfinite(z) for z in (x,y,xm,ym)): continue
                    scales = proj_plane_pixel_scales(w)*3600
                    scale = float(max(abs(scales[0]),abs(scales[1])))
                    cat_inside = 0 <= x < nx and 0 <= y < ny
                    mov_inside = 0 <= xm < nx and 0 <= ym < ny
                    margin_pix = max(1., 1./max(scale,1e-5))
                    near = (-margin_pix<=xm<nx+margin_pix and
                            -margin_pix<=ym<ny+margin_pix)
                    status = ("inside_mosaic_image_wcs" if mov_inside else
                              "near_edge_epoch_uncertain" if near else
                              "outside_mosaic_image_wcs_for_proxy")
                    return dict(wcs_status=status,wcs_hdu=hdu.name,
                                x_catalog=round(float(x),2),y_catalog=round(float(y),2),
                                x_motion_proxy=round(float(xm),2),
                                y_motion_proxy=round(float(ym),2),
                                naxis1=nx,naxis2=ny,pixel_scale_arcsec=round(scale,5),
                                catalog_inbounds=bool(cat_inside),
                                motion_proxy_inbounds=bool(mov_inside),
                                quality_note="Mosaic FITS image WCS only, no detector-level gWCS, SCI/DQ, PSF or detection validation. PM reference year 2023.99 is approximate.",
                                error="")
                except Exception:
                    continue
            return dict(wcs_status="unknown_no_usable_fits_wcs",
                        error="No celestial FITS image WCS; NOT proof of missing target")
    finally:
        signal.alarm(0)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",default="results")
    ap.add_argument("--target",default="target.json")
    ap.add_argument("--limit",type=int,default=18)
    ap.add_argument("--timeout",type=int,default=35)
    args=ap.parse_args()
    root=Path(args.out)
    config=json.loads(Path(args.target).read_text())
    ra,dec=float(config["ra_deg"]),float(config["dec_deg"])
    observations={r["obsid"]:r for r in csv.DictReader((root/"observations.csv").open())
                  if r.get("archive_footprint")=="inside_archive_footprint"}
    products=list(csv.DictReader((root/"science_products.csv").open()))
    selected=[]
    for p in products:
        obs=observations.get(p.get("query_obsid",""))
        # MAST product obsID often differs from the parent archive obsid.
        if not obs: continue
        if not p.get("filename","").lower().endswith(SUFFIXES): continue
        if p.get("rights","").lower() not in ("","public"): continue
        selected.append((priority(obs),obs,p))
    # Audit one representative full visit mosaic per candidate obsid first;
    # then visit-level detector mosaics. This avoids wasting the first
    # batch on dozens of dither/module files from a single visit.
    def audit_rank(entry):\n        priority_group, obs, product = entry\n        name = product.get("filename", "").lower()\n        is_visit_mosaic = name.startswith(obs.get("obs_id", "").lower())\n        return (priority_group, 0 if is_visit_mosaic else 1,\n                obs.get("date_start_utc", ""), name)\n    selected.sort(key=audit_rank)
    unique_first=[]\n    seen_obs=set()\n    for entry in selected:
        key=entry[1].get("obsid", "")\n        if key not in seen_obs:\n            unique_first.append(entry)
            seen_obs.add(key)
    selected=unique_first+[entry for entry in selected if entry not in unique_first]
    path=root/"mosaic_wcs_audit.csv"
    old=list(csv.DictReader(path.open())) if path.exists() else []
    indexed={x["data_uri"]:x for x in old}
    checked=0
    signal.signal(signal.SIGALRM,timeout_handler)
    for _,obs,p in selected:
        uri=p.get("data_uri","")
        if not uri or uri in indexed: continue
        rec=dict(mission=obs.get("mission",""),
                 instrument=obs.get("instrument",""),
                 proposal=obs.get("proposal",""),
                 filter=obs.get("filter",""),
                 obsid=obs.get("obsid",""),
                 obs_id=obs.get("obs_id",""),
                 date_utc=obs.get("date_start_utc",""),
                 product=p.get("filename",""),data_uri=uri)
        try: rec.update(check(uri,ra,dec,obs,args.timeout))
        except Exception as exc:
            rec.update(wcs_status="remote_error_or_timeout",
                       error=(type(exc).__name__+": "+str(exc))[:220])
        indexed[uri]=rec
        checked+=1
        if checked>=args.limit: break
    rows=sorted(indexed.values(),key=lambda r:(priority(r),r.get("date_utc",""),r.get("product","")))
    with path.open("w",newline="") as handle:
        w=csv.DictWriter(handle,fieldnames=FIELDS,extrasaction="ignore")
        w.writeheader();w.writerows(rows)
    counts={}
    for r in rows:counts[r["wcs_status"]]=counts.get(r["wcs_status"],0)+1
    report=dict(updated_utc=datetime.now(timezone.utc).isoformat(),
                audited_mosaic_products=len(rows),checked_this_run=checked,
                available_mosaic_products=len(selected),
                status_counts=counts,source_detections_confirmed=False,
                caveat="Mosaic WCS screening is not source detection, detector-cal gWCS, nor an independent parallax epoch. Approximate PM reference epoch 2023.99.")
    (root/"mosaic_wcs_summary.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))

if __name__=="__main__":main()
