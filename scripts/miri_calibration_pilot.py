#!/usr/bin/env python3
"""MIRI calibrated-detector forced-photometry *diagnostic* for Shabgard.

Each _cal.fits is read with JWST gWCS + SCI/ERR/DQ. MIRI calibrated pixels
are usually MJy/sr; sum of background-subtracted pixel brightnesses is
converted to uJy with the header's pixel solid angle, when available.
Flux is a RAW APERTURE measurement (NO PSF/aperture correction): it must
NOT be mixed with COSMOS2025 AUTO flux or used for silicate inference.

A set of displaced blank apertures supplies a correlated-background proxy.
Does NOT assert identification, detection, publication-ready photometry
or a validated astrometric epoch; no fixed catalog reference MJD exists yet.
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
from audit_calibrated_detector_local import exposure_mjd, sky_at

# Values in reference frame UTC MJD. The COSMOS2025 coordinate epoch is
# UNKNOWN: these are candidate hypotheses, not a posterior/confidence bound.
REF_HYPOTHESES = (60061.0, 60308.0, 60780.0)
AP_RADII = {"F770W": 3.0, "F1000W": 3.5,
            "F1280W": 4.5, "F1500W": 5.0,
            "F1800W": 5.5, "F2100W": 6.0}


def sky_samples(mjd, wcs):
    result = []
    for t0 in REF_HYPOTHESES:
        ra, dec = sky_at(mjd, t0)
        x, y = wcs.invert(ra, dec)
        result.append((t0, float(x), float(y)))
    return result


def local_aperture(sci, err, dq, x, y, radius, save=None):
    ny, nx = sci.shape
    rmax = max(35, math.ceil(radius * 5))
    ix, iy = int(round(x)), int(round(y))
    if not (rmax <= ix < nx-rmax and rmax <= iy < ny-rmax):
        raise ValueError("Position outside detector safe margin for local sky/control test")
    sl_y, sl_x = slice(iy-rmax, iy+rmax+1), slice(ix-rmax, ix+rmax+1)
    a = np.asarray(sci[sl_y, sl_x], dtype=float)
    e = np.asarray(err[sl_y, sl_x], dtype=float)
    q = np.asarray(dq[sl_y, sl_x], dtype=np.uint32)
    yy, xx = np.indices(a.shape, dtype=float)
    xx += ix-rmax
    yy += iy-rmax
    # Only DO_NOT_USE bit (bit 0). Retain all other flags for later audit.
    good = np.isfinite(a) & np.isfinite(e) & (e > 0) & ((q & 1) == 0)
    def measure(cx, cy):
        rr = np.hypot(xx-cx, yy-cy)
        aperture = (rr <= radius) & good
        annulus = (rr >= max(radius * 2.1, radius + 3)) & (rr <= radius * 3.0+4) & good
        npix = int(aperture.sum())
        if npix < max(8, 0.7*math.pi*radius**2) or annulus.sum() < 55:
            return None
        background = float(np.median(a[annulus]))
        signal = float(np.sum(a[aperture]-background))
        formal = float(np.sqrt(np.sum(e[aperture]**2)))
        return signal, formal, background, npix
    found = measure(x, y)
    if found is None:
        raise ValueError("Invalid or missing SCI/DQ/ERR aperture/annulus")
    signal, sigma_err, bkg, n_pix = found
    blank = []
    # Independent-ish spatial samples, avoiding target aperture and annulus.
    # This is a local null estimate, not an injection/recovery experiment.
    for dr in (radius*4.2, radius*5.6, radius*7.0):
        for phase in np.linspace(0, 2*np.pi, 12, endpoint=False):
            u, v = x+dr*np.cos(phase), y+dr*np.sin(phase)
            sample = measure(u, v)
            if sample is not None:
                blank.append(sample[0])
    if len(blank) < 12:
        raise ValueError("Too few valid displaced blank apertures")
    blank = np.asarray(blank)
    blank_offset = float(np.median(blank))
    sig_blank = float(1.4826*np.median(np.abs(blank-blank_offset)))
    # This conservative proxy does not remove covariance within ERR pixels.
    sigma = max(sigma_err, sig_blank)
    if sigma <= 0 or not np.isfinite(sigma):
        raise ValueError("Invalid photometric uncertainty proxy")
    if save is not None:
        save.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(save, sci=a.astype(np.float32), err=e.astype(np.float32),
                            dq=q, x_target=x-(ix-rmax), y_target=y-(iy-rmax),
                            aperture_radius=radius, background=bkg,
                            null_aperture_signals=blank)
    return dict(aperture_raw_sum_mjysr=signal, sigma_err_sum_mjysr=sigma_err,
                sigma_blank_mjysr=sig_blank, sigma_proxy_mjysr=sigma,
                aperture_signal_to_noise_proxy=(signal-blank_offset)/sigma,
                aperture_n_good_pix=n_pix, local_background_mjysr=bkg,
                n_blank_controls=len(blank), blank_aperture_median_mjysr=blank_offset,
                blank_controls_ge_target=int(np.sum(blank >= signal)))


def inspect(path, outdir):
    row = dict(filename=path.name, status="error", reason="",
               source_detected=False, astrometric_epoch_usable=False,
               photometry_publication_ready=False)
    try:
        if not path.name.endswith("_mirimage_cal.fits"):
            raise ValueError("Not a MIRI detector calibrated _cal.fits file")
        from jwst import datamodels
        with fits.open(path, memmap=False) as f, datamodels.open(str(path)) as model:
            sci, err, dq = f["SCI"].data, f["ERR"].data, f["DQ"].data
            p = f[0].header
            s = f["SCI"].header
            filt = str(p.get("FILTER", model.meta.instrument.filter))
            mjd = exposure_mjd(p, s)
            wcs = model.meta.wcs
            if wcs is None:
                raise ValueError("Missing JWST detector gWCS")
            row.update(filter=filt, mid_mjd=mjd, sci_bunit=s.get("BUNIT", ""),
                       detector=p.get("DETECTOR", ""), instrument=p.get("INSTRUME", ""))
            refs = sky_samples(mjd, wcs)
            for j, (_, x, y) in enumerate(refs):
                row["ref_mjd_"+str(j)] = refs[j][0]
                row["x_ref_"+str(j)] = x
                row["y_ref_"+str(j)] = y
            ny, nx = sci.shape
            rad = AP_RADII.get(filt, 4)
            inside = [np.isfinite(x+y) and rad+2 <= x < nx-rad-2
                      and rad+2 <= y < ny-rad-2 for _, x, y in refs]
            row["n_ref_hypotheses_inside"] = sum(inside)
            row["coverage_any_ref_hypothesis"] = any(inside)
            row["coverage_all_ref_hypotheses"] = all(inside)
            if not any(inside):
                row["status"] = "outside_all_tested_reference_positions"
                row["reason"] = ("Only 3 provisional catalog epochs tested; "
                                 "no formal absence-of-coverage claim")
                return row
            _, x, y = refs[1] if inside[1] else refs[inside.index(True)]
            if not inside[1]:
                row["reason"] = "Nominal reference epoch outside: testing alternative only"
            row.update(x_photometry=x, y_photometry=y, aperture_radius_pix=rad)
            try:
                vals = local_aperture(sci, err, dq, x, y, rad,
                                      outdir/"cutouts"/(path.stem+"_local.npz"))
            except ValueError as e:
                row["status"] = "inside_gwcs_but_photometry_unusable"
                row["reason"] = str(e)
                return row
            row.update(vals)
            bunit = str(s.get("BUNIT", "")).lower().replace(" ", "")
            area = p.get("PIXAR_SR", s.get("PIXAR_SR", None))
            if area is None:
                area = getattr(model.meta.photometry, "pixelarea_steradians", None)
            row["pixel_area_sr"] = str(area) if area is not None else ""
            if bunit in ("mjy/sr", "mjy/sr-1") and area is not None:
                k = float(area) * 1e12
                row["flux_aperture_uncorrected_ujy"] = vals["aperture_raw_sum_mjysr"] * k
                row["flux_sigma_proxy_ujy"] = vals["sigma_proxy_mjysr"] * k
            else:
                row["reason"] = "Missing or unrecognized BUNIT/PIXAR_SR; native-unit proxy only"
            row["status"] = "detector_pixel_and_forced_aperture_DIAGNOSTIC"
            row["interpretation"] = (
                "NOT source identification: no MIRI PSF/aperture correction, "
                "neighbor deblending or image-to-image stationary registration; "
                "blank controls are NOT an injection-recovery test")
            return row
    except Exception as e:
        row["reason"] = type(e).__name__ + ": " + str(e)[:260]
        return row


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("files", nargs="+", help="Glob patterns for MIRI CAL files")
    p.add_argument("--out", default="miri_results")
    args = p.parse_args()
    paths = [Path(x) for pattern in args.files for x in sorted(glob.glob(pattern))]
    if not paths:
        p.error("No _mirimage_cal.fits inputs matched")
    output = Path(args.out)
    output.mkdir(parents=True, exist_ok=True)
    rows = [inspect(path, output) for path in paths]
    fields = sorted(set().union(*(r.keys() for r in rows)))
    with (output/"miri_cal_detector_diagnostic.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    from collections import Counter
    summary = dict(n_files=len(rows),
                   status_counts=dict(Counter(r["status"] for r in rows)),
                   filter_counts=dict(Counter(r.get("filter", "unknown") for r in rows)),
                   n_with_gwcs_coverage=sum(bool(r.get("coverage_any_ref_hypothesis")) for r in rows),
                   n_diagnostic_aperture=sum(r["status"] == "detector_pixel_and_forced_aperture_DIAGNOSTIC" for r in rows),
                   source_detected=False, silicate_feature_detected=False,
                   astrometric_epochs_validated=0,
                   caveat="Provisional coordinate epochs, no PSF or aperture correction; cannot claim SED/colors from diagnostic CAL aperture sums")
    (output/"miri_cal_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))
    # Fail the workflow when code/inputs fail everywhere, so Actions cannot
    # misleadingly report scientific success after a full analysis error.
    if all(r["status"] == "error" for r in rows):
        raise SystemExit("All MIRI detector checks errored; inspect CSV and logs")


if __name__ == "__main__":
    main()
