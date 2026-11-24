#!/usr/bin/env python3
"""Robust Gaussian localization on MIRI CAL cutout pixels (DIAGNOSTIC ONLY).

Useful for source association and inter-dither consistency checks.
A MIRI PSF is NOT a Gaussian. No stationary-source registration or
COSMOS2025 effective coordinate epoch has yet been verified; do not
treat formal uncertainties here as absolute astrometric accuracy.
"""
import argparse
import csv
import glob
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


def inspect(path, radius=5):
    r = {"cutout": str(path), "status": "error", "reason": ""}
    try:
        with np.load(path) as d:
            a = d["sci"].astype(float)
            e = d["err"].astype(float)
            dq = d["dq"].astype(np.uint32)
            x = float(d["x_target"])
            y = float(d["y_target"])
            bg = float(d["background"])
        yy, xx = np.indices(a.shape, dtype=float)
        d2 = np.hypot(xx-x, yy-y)
        # Conservative for source fit: discard ALL flagged pixels.
        good = np.isfinite(a) & np.isfinite(e) & (e > 0) & (dq == 0) & (d2 <= radius)
        if good.sum() < 50:
            raise ValueError("Fewer than 50 clean local pixels")
        X, Y, flux = (xx[good]-x), (yy[good]-y), a[good]
        err = np.maximum(e[good], np.median(e[good]))
        init = [max(float(np.percentile(flux, 90)-bg), .01), 0, 0, 1.0, bg, 0, 0]
        def model(p):
            amp, dx, dy, sigma, b, bx, by = p
            return (amp * np.exp(-.5*((X-dx)**2+(Y-dy)**2)/sigma**2)
                    + b + bx*X + by*Y)
        low = [0, -3, -3, .6, -np.inf, -np.inf, -np.inf]
        high = [np.inf, 3, 3, 2.6, np.inf, np.inf, np.inf]
        solution = least_squares(lambda p: (model(p)-flux)/err, init,
                                 bounds=(low, high), loss="soft_l1",
                                 f_scale=2, max_nfev=500)
        if not solution.success:
            raise ValueError("Fit failed: "+solution.message)
        p = solution.x
        redchi = float(np.sum(((model(p)-flux)/err)**2) / (len(flux)-7))
        formal = np.full(7, np.nan)
        if np.linalg.cond(solution.jac.T @ solution.jac) < 1e12:
            cov = np.linalg.pinv(solution.jac.T @ solution.jac) * max(redchi, 1.)
            formal = np.sqrt(np.maximum(np.diag(cov), 0))
        r.update(status="gaussian_candidate_diagnostic", x_offset_pix=float(p[1]),
                 y_offset_pix=float(p[2]),
                 x_error_formal_pix=float(formal[1]), y_error_formal_pix=float(formal[2]),
                 sigma_gaussian_pix=float(p[3]), amplitude_native=float(p[0]),
                 amplitude_over_formal_error=float(p[0]/formal[0]) if formal[0]>0 else "",
                 reduced_chi2=redchi, pixels_fit=int(good.sum()),
                 at_centroid_bound=bool(abs(p[1])>2.98 or abs(p[2])>2.98))
    except Exception as exc:
        r["reason"] = type(exc).__name__+": "+str(exc)[:220]
    return r


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("cutouts", nargs="+")
    ap.add_argument("--out", default="miri_results/miri_gaussian_centroids.csv")
    args=ap.parse_args()
    paths = [Path(z) for pat in args.cutouts for z in sorted(glob.glob(pat))]
    if not paths:
        raise SystemExit("No cutouts matched; only CAL-gWCS-positive cutouts are analyzed")
    rows = [inspect(p) for p in paths]
    dest = Path(args.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    cols = sorted({k for r in rows for k in r.keys()})
    with dest.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=cols)
        writer.writeheader()
        writer.writerows(rows)
    ok = [r for r in rows if r["status"]=="gaussian_candidate_diagnostic"]
    summary = {"n_cutouts":len(rows), "n_fits":len(ok),
               "reference_epoch_verified":False,
               "astrometric_detection_confirmed":False,
               "miri_psf_fit_completed":False,
               "method":"Gaussian + robust plane background, diagnostic only",
               "caveat":"No ePSF, stationary-field registration, or independently verified position provenance"}
    if len(ok)>1:
        summary["centroid_x_offset_sample_sd_pix"]=float(np.std([x["x_offset_pix"] for x in ok],ddof=1))
        summary["centroid_y_offset_sample_sd_pix"]=float(np.std([x["y_offset_pix"] for x in ok],ddof=1))
    (dest.parent/"miri_gaussian_centroids_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))
    if len(ok)!=len(rows):
        raise SystemExit("Some MIRI centroid diagnostic fits failed; inspect CSV")


if __name__=="__main__":
    main()
