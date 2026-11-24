# Shabgard MIRI calibrated-detector science checkpoint (2026-11-24)

Target: COSMOS2025 ID 639188, RA 150.0634341, Dec +2.3552216.
Approximate PM (mu_alpha*cos(delta), mu_delta)=(-70,-54) mas/yr.
**COSMOS2025 effective catalog coordinate reference epoch remains unverified.**

## Reproducible input evidence

- [MIRI pilot run 45410470478](https://github.com/physicsbahare/Shabgard/actions/runs/45410470478), artifact `shabgard-miri-detector-pilot`: initial F770W CAL coverage and diagnostic aperture measurements.
- [All-pointing scan 45411234868](https://github.com/physicsbahare/Shabgard/actions/runs/45411234868), artifact `shabgard-miri-pointing-coverage`: F1000W, F1280W, F1500W CAL pointing identification. **Not all dithers were inspected in this scan.**
- [Covered-dither followup run 45411830409](https://github.com/physicsbahare/Shabgard/actions/runs/45411830409), artifact `shabgard-miri-covered-dithers`: 23 real MIRI CAL files, 23 detector-gWCS geometric positions inside, 21 usable SCI/ERR/DQ local aperture diagnostics and 21 robust-Gaussian candidate fits.

All quantities below are from `miri_cal_detector_diagnostic.csv` and `miri_gaussian_centroids.csv` in the **covered-dither followup artifact**, not archive-footprint estimates.

| Filter | UTC date | MIRI CAL exposures with local usable SCI/ERR/DQ | Uncorrected aperture-flux median (uJy) | SD across exposures (uJy) | Median single-exposure local sigma proxy (uJy) | Gaussian x/y centroid-offset SD relative to a PROVISIONAL predicted position (pixels) |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| F770W | 2024-01-04 | 8 | 0.6064 | 0.0968 | 0.0948 | 0.125 / 0.119 |
| F1000W | 2024-12-21 | 1 of 3 gWCS-positive exposures (2 invalid SCI/ERR/DQ annulus/aperture) | 1.5327 | not estimable | 0.1990 | not estimable |
| F1280W | 2025-12-29 | 6 | 1.6162 | 0.4105 | 0.3225 | 0.524 / 1.010 |
| F1500W | 2025-12-29 | 6 | 1.6119 | 0.2765 | 0.5331 | 0.547 / 0.940 |

**Interpretation restrictions:**
- These are **uncorrected aperture fluxes**, not total PSF-corrected fluxes, synthetic-band photometry, or published detections. Apertures differ with wavelength; comparing their raw flux ratios **cannot** test the 10 μm silicate feature.
- Background estimates use local annuli and displaced blank apertures. Their noise estimates are proxies; residual nearby sources, correlated detector noise and PSF wings need explicit treatment. Only one F1000W exposure passed local SCI/ERR/DQ checks.
- Robust circular-Gaussian fits are diagnostic, NOT MIRI ePSF/physical-PSF astrometry. Position offsets use provisional COSMOS2025 coordinate reference MJD 60308. No stationary-galaxy registration, joint source-association or independent parallax factors have been performed.
- F770W Gaussian fit offsets have mean approximately (-1.367,-0.497) detector pixels vs provisional position, with inter-exposure SD (0.125,0.119) pixels. This inter-dither consistency is promising, but is not absolute RA/Dec accuracy or proof of identity.
- The COSMOS2025 catalog contains `snr_f770w=143.1589858`, `flux_auto_f770w=0.95145434`, but catalog units/source matching and PSF-to-AUTO aperture corrections **must** be checked before quantitative comparison.
- These visits are distinct MIRI observing dates, **not** yet validated astrometric epochs. No 5-parameter parallax fit or silicate identification is scientifically justified from this pipeline.

## Next needed science (ordered)

1. Validate candidate source association in F770W with measured PM and fixed stationary field control sources; determine COSMOS2025 effective coordinate epoch from B7 detection mosaic contributors.
2. Fit physical or empirical MIRI PSFs plus nearby-source model in all 4 bands, using robust background and calibrated pixel covariance; inject/recover simulated sources and test blank positions. Obtain aperture-corrected fluxes/upper limits with uncertainty budgets and verify `BUNIT`/`PIXAR_SR` and F770W catalog matching.
3. For F1000W two gWCS-positive-but-unusable exposures, inspect local SCI/ERR/DQ masks individually; do not assume a nondetection. Inspect remaining candidate CAL pointings/exposures only as justified by footprint.
4. Convolve suitable ~720 K metal-poor cloudy and cloud-free substellar SEDs with **real STScI MIRI throughput curves**; test 10 μm silicate vs NH3/CH4/H2O alternatives with time variability/systematics. Photometry cannot resolve mineralogical detail like MIRI spectroscopy.
5. Separately quantify astrometric reference registration and true spacecraft parallax factors before considering any MIRI visit as a parallax epoch.

Source scripts: `scripts/select_miri_cal_pilot.py`, `scripts/download_miri_cal_batch.py`, `scripts/miri_calibration_pilot.py`, `scripts/miri_cutout_centroids.py`; workflows `.github/workflows/miri-cal-pilot.yml`, `.github/workflows/miri-cal-coverage-scan.yml`, `.github/workflows/miri-cal-followup.yml`.
