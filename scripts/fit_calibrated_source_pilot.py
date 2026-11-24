#!/usr/bin/env python3
"""Pilot centroids on JWST NIRCam calibrated exposures, NOT publishable astrometry.

Fits a small, deliberately approximate elliptical-Gaussian source plus planar
background to each independent detector SCI/ERR/DQ cutout. Real JWST PSFs are
NOT Gaussian: these fits only test candidate-source positional consistency
before empirical/ePSF, stationary-reference registration, and parallax fitting.

The supplied COSMOS2025 coordinate reference MJD is a hypothesis until
the detection-coadd provenance has been verified.
"""
import argparse
import csv
import glob
import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.optimize import least_squares
from astropy.time import Time

RA, DEC = 150.0634341, 2.3552216
MU_RA, MU_DEC = -70.0, -54.0  # mas/year in (alpha*cos(dec), delta)


def midpoint(primary, sci):
    for h in (sci, primary):
        for k in ("MJD-AVG", "MJD_AVG", "MJD-MID", "MJD_MID"):
            if k in h:
                return float(h[k])
    for h in (sci, primary):
        if "EXPSTART" in h and "EXPEND" in h:
            return (float(h["EXPSTART"]) + float(h["EXPEND"])) / 2.0
    for h in (sci, primary):
        if "EXPSTART" in h and "EXPTIME" in h:
            return float(h["EXPSTART"]) + float(h["EXPTIME"]) / 172800.0
    for h in (sci, primary):
        if "DATE-OBS" in h and "EXPTIME" in h:
            st = h["DATE-OBS"] + ("T" + h["TIME-OBS"] if "TIME-OBS" in h else "")
            return Time(st).mjd + float(h["EXPTIME"]) / 172800.0
    raise ValueError("Exposure midpoint cannot be determined from FITS headers")


def predict(mjd, ref_mjd):
    dt = (mjd - ref_mjd) / 365.25
    return (RA + MU_RA * dt / (3.6e6 * math.cos(math.radians(DEC))),
            DEC + MU_DEC * dt / 3.6e6)


def gaussian_plane(p, xx, yy, centerx, centery):
    amp, dx, dy, sx, sy, background, bgx, bgy = p
    xp = (xx - centerx) - dx
    yp = (yy - centery) - dy
    return (amp * np.exp(-0.5 * ((xp/sx)**2 + (yp/sy)**2))
            + background + bgx * (xx - centerx) + bgy * (yy - centery))


def fit_patch(data, error, dq, px, py, radius=9):
    ny, nx = data.shape
    # all sampled pixels must lie on the detector and have real SCI/ERR/DQ
    ix, iy = round(px), round(py)
    if (ix < radius+1 or iy < radius+1
            or ix >= nx-radius-1 or iy >= ny-radius-1):
        raise ValueError("Target predicted position outside detector safe margin")
    cut = np.asarray(data[iy-radius:iy+radius+1, ix-radius:ix+radius+1],
                     dtype=np.float64)
    if error is not None:
        err = np.asarray(error[iy-radius:iy+radius+1, ix-radius:ix+radius+1],
                         dtype=np.float64)
    else:
        err = None
    yy, xx = np.indices(cut.shape, dtype=float)
    xx += ix-radius
    yy += iy-radius
    mask = np.isfinite(cut)
    if dq is not None:
        mask &= (np.asarray(dq[iy-radius:iy+radius+1,
                               ix-radius:ix+radius+1]) == 0)
    if err is not None:
        mask &= np.isfinite(err) & (err > 0)
    annulus = (np.hypot(xx-px, yy-py) >= 5.5) & mask
    if mask.sum() < 120 or annulus.sum() < 40:
        raise ValueError("Too few good SCI/DQ/ERR or local sky pixels")
    bkg = float(np.median(cut[annulus]))
    local_sigma = 1.4826 * float(np.median(np.abs(cut[annulus]-bkg)))
    if not np.isfinite(local_sigma) or local_sigma <= 0:
        raise ValueError("Nonpositive local sky MAD")
    if err is None:
        weights = np.full_like(cut, local_sigma)
        noise_source = "sky_MAD_only"
    else:
        # ERR is the pipeline per-pixel random uncertainty; it is not the
        # registration or PSF-model systematic. Fit a floor from sky MAD.
        weights = np.maximum(err, local_sigma * 0.3)
        noise_source = "ERR_plus_MAD_floor"
    # Fit at an independently predicted position; do not allow a centroid
    # optimizer to latch onto an unrelated source far from that position.
    amp0 = max(float(np.nanmax(cut[mask]) - bkg), local_sigma)
    guess = np.array([amp0, 0., 0., 0.9, 0.9, bkg, 0., 0.])
    lower = np.array([0., -2.5, -2.5, .55, .55, -np.inf, -np.inf, -np.inf])
    upper = np.array([np.inf, 2.5, 2.5, 2.2, 2.2, np.inf, np.inf, np.inf])
    def residual(par):
        return ((gaussian_plane(par, xx, yy, px, py) - cut) / weights)[mask]
    solution = least_squares(residual, guess, bounds=(lower, upper),
                             loss="soft_l1", f_scale=2., max_nfev=250)
    if not solution.success:
        raise ValueError("Least squares did not converge: " + solution.message)
    par = solution.x
    dof = int(mask.sum() - len(par))
    r = residual(par)
    chi2_reduced = float(np.sum(r*r) / dof) if dof > 0 else float("nan")
    # Covariance is FORMAL ONLY; soft-loss, noise correlations, PSF mismatch
    # and positional calibration mean it will underestimate true uncertainties.
    cx_err = cy_err = amp_err = float("nan")
    try:
        jtj = solution.jac.T @ solution.jac
        if np.linalg.cond(jtj) < 1e12:
            cov = np.linalg.inv(jtj) * max(1., chi2_reduced)
            if np.all(np.diag(cov) >= 0):
                amp_err, cx_err, cy_err = np.sqrt(np.diag(cov))[[0,1,2]]
    except np.linalg.LinAlgError:
        pass
    on_bounds = any(abs(par[i] - lo) < 0.005 or abs(hi-par[i]) < 0.005
                    for i, (lo, hi) in enumerate(zip(lower, upper)) if
                    np.isfinite(lo) and np.isfinite(hi))
    # Derived fit scores are NOT independent source-detection significance.
    score = float(par[0]/amp_err) if amp_err > 0 and np.isfinite(amp_err) else float("nan")
    return {
        "x_fit": float(px+par[1]), "y_fit": float(py+par[2]),
        "x_err_formal_pix": float(cx_err), "y_err_formal_pix": float(cy_err),
        "amplitude": float(par[0]), "amplitude_err_formal": float(amp_err),
        "amplitude_over_formal_error": score, "sigma_x_pix": float(par[3]),
        "sigma_y_pix": float(par[4]), "background": float(par[5]),
        "local_sky_mad": float(local_sigma), "reduced_chi2": chi2_reduced,
        "centroid_offset_pix": float(np.hypot(par[1], par[2])),
        "fit_at_parameter_limit": bool(on_bounds),
        "good_pixels": int(mask.sum()), "noise_source": noise_source,
        "method": "elliptical_Gaussian_DIAGNOSTIC_only",
    }


def inspect(path, reference_mjd):
    out = dict(filename=path.name, status="error", reason="")
    try:
        if not path.name.lower().endswith("_cal.fits"):
            raise ValueError("Not a JWST calibrated detector _cal.fits")
        with fits.open(path, memmap=False) as h:
            p = h[0].header
            sci = h["SCI"]
            mjd = midpoint(p, sci.header)
            ra, dec = predict(mjd, reference_mjd)
            out.update(mid_mjd=mjd, ra_pred=ra, dec_pred=dec,
                       filter=str(p.get("FILTER", sci.header.get("FILTER", "UNKNOWN"))))
            err = h["ERR"].data if "ERR" in h else None
            dq = h["DQ"].data if "DQ" in h else None
            if err is None or dq is None:
                raise ValueError("Missing ERR or DQ: cannot use as photometric pilot")
            # JWST detector-level FITS uses gWCS stored in ASDF, not a generic
            # linear FITS WCS. Always use JWST datamodels for this transformation.
            from jwst import datamodels
            with datamodels.open(str(path)) as model:
                if model.meta.wcs is None:
                    raise ValueError("JWST gWCS absent")
                px, py = model.meta.wcs.invert(ra, dec)
                px, py = float(px), float(py)
                if not np.isfinite(px+py):
                    raise ValueError("Non-finite sky-to-detector gWCS")
                out.update(x_pred=px, y_pred=py)
                fit = fit_patch(sci.data, err, dq, px, py)
                out.update(fit)
                radec = model.meta.wcs(fit["x_fit"], fit["y_fit"])
                fit_ra, fit_dec = float(radec[0]), float(radec[1])
                out.update(ra_fit=fit_ra, dec_fit=fit_dec,
                           delta_ra_cosdec_mas=(fit_ra-ra)*3.6e6*math.cos(math.radians(dec)),
                           delta_dec_mas=(fit_dec-dec)*3.6e6)
                out.update(status="fit_diagnostic_only",
                           reason="Gaussian is NOT a JWST PSF; no stationary-control "
                                  "registration, field ePSF, catalog epoch or parallax factors")
    except Exception as e:
        out["reason"] = type(e).__name__ + ": " + str(e)[:240]
    return out


def clean(x):
    if isinstance(x, (np.floating, float)) and not np.isfinite(x):
        return None
    if isinstance(x, np.integer):
        return int(x)
    return x


def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument("files", nargs="+")
    a.add_argument("--reference-mjd", type=float, required=True,
                   help="PROVISIONAL catalog coordinate reference hypothesis")
    a.add_argument("--out", default="pilot_results/gaussian_centroid_diagnostic.csv")
    args = a.parse_args()
    paths = [Path(p) for pattern in args.files for p in sorted(glob.glob(pattern))]
    if not paths:
        a.error("No calibrated FITS matches")
    rows = [inspect(p, args.reference_mjd) for p in paths]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    columns = sorted({k for row in rows for k in row})
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    good = [r for r in rows if r["status"] == "fit_diagnostic_only"]
    summary = {
        "n_files": len(rows), "n_diagnostic_fits": len(good),
        "catalog_ref_mjd_assumed": args.reference_mjd,
        "reference_epoch_verified": False,
        "real_detection_confirmed": False, "usable_astrometric_epoch_confirmed": False,
        "reason": "Gaussian source-localization pilot only; JWST PSF, reference stars/galaxies, "
                  "astrometric registration and parallax factors are missing.",
        "dither_ra_residual_rms_mas": None,
        "dither_dec_residual_rms_mas": None,
    }
    if len(good) > 1:
        ra = np.asarray([r["delta_ra_cosdec_mas"] for r in good])
        dec = np.asarray([r["delta_dec_mas"] for r in good])
        summary["dither_ra_residual_rms_mas"] = float(np.std(ra, ddof=1))
        summary["dither_dec_residual_rms_mas"] = float(np.std(dec, ddof=1))
    (out.parent / "gaussian_centroid_diagnostic_summary.json").write_text(
        json.dumps(summary, indent=2, default=clean) + "\n")
    print(json.dumps(summary, indent=2, default=clean))
    if len(good) != len(rows):
        raise SystemExit("Some calibrated detector fits failed; inspect output artifacts")


if __name__ == "__main__":
    main()
