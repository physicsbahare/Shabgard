# Shabgard COSMOS2025 coordinate epoch: evidence and working policy

Updated 2026-11-21. Target COSMOS2025 ID 639188, RA=150.06343413982992 deg, Dec=+2.355221569555261 deg. Approximate published proper motion (mu_alpha*cos(dec), mu_dec)=(-70,-54) mas/yr. **Do not describe an assumed epoch as measured.**

## Evidence actually present in this repository

- `data/shabgard_cosmos2025_catalog_info.txt` contains the complete extracted target row and eight FITS headers. It explicitly reports **no coordinate reference epoch** for the catalog RA/Dec columns.
- `results/epoch_groups.csv` lists JWST PID 1727, observation `jw01727-o051_t051_nircam_clear-f115w`, `...f150w`, `...f277w`, and `...f444w` on **2023-12-30**, each as an archive-footprint candidate. These four entries are not detector-level proof of coverage, nor proof that the COSMOS2025 detection image or measured catalog coordinates correspond to exactly that date.
- The same file lists PID 1727 MIRI F770W archive-footprint candidates on 2024-01-04 (o016) and 2024-01-06 (o018); these are **not** the four-band NIRCam detection reference.
- COSMOS2025 tile B7 is associated with the source in the target extraction, but an explicit detector/exposure-to-tile-B7 provenance and effective detection-image MJD remain to be verified.

## Operational reference date (provisional, NOT catalog metadata)

Use **2023-12-30** (UTC date only; exact exposure mid-MJD still unknown) as a *working hypothesis* for the COSMOS-Web NIRCam coordinate measurement. This is stronger than guessing from the 2025 catalog release date, but **not a validated catalog astrometric reference epoch**. Do not hard-code 2023.99 or 2023-12-30 as a confirmed epoch. If source positions were measured on a multi-visit coadd, there may be no unique exposure epoch; the effective epoch depends on contributing visits and their weights.

For approximate position propagation, with t and t0 in Julian years:
- delta(alpha*cos(delta)) = -70*(t-t0) mas
- delta(delta) = -54*(t-t0) mas
- alpha(t) = alpha0 + delta(alpha*cos(delta))/(3.6e6*cos(delta0)) deg
- delta(t) = delta0 + delta(delta)/3.6e6 deg

For a **2011-12-10** HST exposure and provisional t0~2023.9945, dt~-12.06 yr, expected position offset relative to catalog is **+0.844 arcsec east** and **+0.651 arcsec north**, total ~1.066 arcsec, *before parallax*. These numbers are approximate, not measured HST detections. For a one-year reference-epoch error, the position shifts ~0.088 arcsec. Do not claim HST noncoverage based on a single unverified propagated point. Use a search region covering plausible reference dates, proper-motion uncertainties, astrometric registration and parallax.

## Required validation, in order

1. Identify the actual PID 1727 o051 calibrated `*_cal.fits` exposures and exposure-level start/end/mid MJD, detectors, filters, and footprints; test target location against each calibrated detector's gWCS and SCI/DQ pixels.
2. Identify exactly which NIRCam exposures/visits feed the COSMOS2025 B7 four-band detection coadd. If multiple dates contribute, calculate an exposure/weight-aware effective date, and document that it is a proxy rather than an official epoch unless authors define one.
3. For each HST WFC3/IR and MIRI candidate, test calibrated detector-level WCS and actual usable SCI/DQ coverage at a proper-motion- and parallax-aware search region; report source detection/SNR and centroid errors separately.
4. For any actual detections, report filter, instrument, exposure mid-MJD, detector, measured RA/Dec, astrometric covariance, and observer parallax factors. Only count independent usable astrometric epochs, not MAST footprints or mosaic image-WCS positives.
5. Do not claim feasibility of a five-parameter (position, two proper motions, parallax) fit until sufficient independently timed, calibrated detections and parallax-factor diversity exist.

## Status

**Confirmed:** a 2023-12-30 NIRCam *archive-footprint candidate*, not the coordinate epoch. **Unresolved:** effective COSMOS2025 B7 coordinate epoch, calibrated detector coverage and actual detections. This is a provenance question, not a reason to stop HST image-level searching.
