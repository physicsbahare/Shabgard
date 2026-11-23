#!/usr/bin/env python3
"""Conservative detector-level Shabgard audit of LOCAL calibrated FITS files.

Requires numpy, astropy; JWST *_cal.fits additionally requires jwst.
Never equates a MAST footprint, i2d mosaic, or pixel availability with detection.
"""
import argparse
import csv
import glob
import math
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u

RA, DEC = 150.0634341, 2.3552216
MU_RA, MU_DEC = -70., -54.  # mas/yr; mu_alpha*cos(delta)


def sky_at(mjd, reference_mjd):
    if reference_mjd is None:
        raise ValueError("Reference MJD is mandatory: COSMOS2025 coordinate epoch unverified")
    dt = (mjd - reference_mjd) / 365.25
    return (RA + MU_RA * dt / (3.6e6 * math.cos(math.radians(DEC))),
            DEC + MU_DEC * dt / 3.6e6)


def exposure_mjd(primary, science):
    for hdr in (science, primary):
        for key in ("MJD-AVG", "MJD_AVG", "MJD-MID", "MJD_MID"):
            if key in hdr:
                return float(hdr[key])
    for hdr in (science, primary):
        if "EXPSTART" in hdr and "EXPEND" in hdr:
            return (float(hdr["EXPSTART"]) + float(hdr["EXPEND"])) / 2.0
    for hdr in (science, primary):
        if "EXPSTART" in hdr and "EXPTIME" in hdr:
            return float(hdr["EXPSTART"]) + float(hdr["EXPTIME"]) / 172800.0
    for hdr in (science, primary):
        if "EXPSTART" in hdr:
            raise ValueError("Only exposure start MJD found; midpoint cannot be verified")
    for hdr in (science, primary):
        date = hdr.get("DATE-OBS")
        if date:
            return float(Time(date + ("T" + hdr["TIME-OBS"] if "TIME-OBS" in hdr else "")).mjd)
    raise ValueError("No exposure time found")


def inspect(path, reference_mjd, radius):
    out = dict(path=str(path), product_type="unknown", mjd="", ra_pred="", dec_pred="",
               pixel_x="", pixel_y="", detector_covered=False, good_pixels=False,
               local_snr="", centroid_x="", centroid_y="", detected=False,
               usable_astrometric_epoch=False, reason="")
    try:
        with fits.open(path, memmap=True) as hdul:
            p = hdul[0].header
            sci_hdus = [h for h in hdul if h.name == "SCI" and h.data is not None]
            if not sci_hdus:
                out["reason"] = "No SCI extension (not a validated calibrated exposure)"
                return out
            sci = sci_hdus[0]
            arr = np.asarray(sci.data)
            if arr.ndim != 2:
                out["reason"] = "SCI is not 2-D"
                return out
            name = path.name.lower()
            jwst_cal = name.endswith("_cal.fits")
            hst_flt = name.endswith(("_flt.fits", "_flc.fits"))
            if not (jwst_cal or hst_flt):
                out["reason"] = "Not a supported individual calibrated cal/flt/flc product"
                return out
            out["product_type"] = "JWST_CAL" if jwst_cal else "HST_FLT_FLC"
            mjd = exposure_mjd(p, sci.header)
            out["mjd"] = round(mjd, 7)
            ra, dec = sky_at(mjd, reference_mjd)
            out["ra_pred"], out["dec_pred"] = ra, dec
            if jwst_cal:
                try:
                    from jwst import datamodels
                except ImportError:
                    out["reason"] = "jwst package required for detector gWCS"
                    return out
                with datamodels.open(str(path)) as model:
                    if model.meta.wcs is None:
                        out["reason"] = "Missing JWST gWCS"
                        return out
                    x, y = model.meta.wcs.invert(ra, dec)
                    x, y = float(x), float(y)
            else:
                wcs = WCS(sci.header)
                x, y = [float(v) for v in wcs.all_world2pix([[ra, dec]], 0)[0]]
            out["pixel_x"], out["pixel_y"] = round(x, 4), round(y, 4)
            ny, nx = arr.shape
            if not (radius + 1 <= x < nx - radius - 1 and radius + 1 <= y < ny - radius - 1):
                out["reason"] = "Outside usable detector SCI margin"
                return out
            out["detector_covered"] = True
            ix, iy = int(round(x)), int(round(y))
            patch = np.asarray(arr[iy-radius:iy+radius+1, ix-radius:ix+radius+1], dtype=float)
            yy, xx = np.indices(patch.shape)
            rr = np.hypot(xx + ix-radius-x, yy + iy-radius-y)
            dq = next((h.data for h in hdul if h.name == "DQ" and h.data is not None), None)
            good = np.isfinite(patch)
            if dq is not None:
                dq_patch = np.asarray(dq[iy-radius:iy+radius+1, ix-radius:ix+radius+1])
                good &= (dq_patch == 0)  # conservative; instrument-specific bit handling later
            inner = (rr <= 2) & good
            annulus = (rr >= 4) & (rr <= radius) & good
            if inner.sum() < 5 or annulus.sum() < 20:
                out["reason"] = "Insufficient valid SCI/DQ pixels for aperture/background"
                return out
            out["good_pixels"] = True
            bkg = float(np.median(patch[annulus]))
            mad = float(np.median(np.abs(patch[annulus] - bkg)))
            sigma = 1.4826 * mad
            if not np.isfinite(sigma) or sigma <= 0:
                out["reason"] = "Invalid local noise estimate"
                return out
            signal = float(np.sum(patch[inner] - bkg))
            snr = signal / (sigma * math.sqrt(int(inner.sum())))
            out["local_snr"] = round(snr, 3)
            out["reason"] = ("Aperture proxy only; no PSF fit, centroid covariance, "
                             "stationary controls, or independent parallax factors")
            # Deliberately do NOT mark detected or usable_astrometric_epoch.
            return out
    except Exception as exc:
        out["reason"] = type(exc).__name__ + ": " + str(exc)[:180]
        return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("files", nargs="+", help="Individual calibrated FITS files (not mosaics)")
    ap.add_argument("--reference-mjd", type=float, required=True,
                    help="Verified epoch of supplied COSMOS2025 RA/Dec, not a guessed value")
    ap.add_argument("--out", default="results/detector_local_audit.csv")
    ap.add_argument("--radius", type=int, default=7)
    args = ap.parse_args()
    if args.radius < 5:
        ap.error("--radius must be at least 5")
    paths = [Path(p) for pattern in args.files for p in sorted(glob.glob(pattern))]
    if not paths:
        ap.error("No FITS files matched")
    rows = [inspect(p, args.reference_mjd, args.radius) for p in paths]
    dest = Path(args.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} calibrated-product checks to {dest}")
    print("WARNING: local_snr is a proxy, NOT a confirmed detection or astrometric epoch")


if __name__ == "__main__":
    main()
