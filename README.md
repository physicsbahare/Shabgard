# Shabgard: all-epoch space-telescope archive census

Goal: **find every public HST and JWST observation that may contain Shabgard**, including NIRCam, MIRI, NIRISS, NIRSpec and HST ACS/WFC3, without mistaking multiple dithers/products for separate parallax epochs.

Target: COSMOS2025 ID **639188**, Shabgard (Salmasi et al. 2026, PASP; arXiv:2610.07861).

## Status (2026-11-16)
Automated MAST discovery is being established. **No new coverage or MIRI detection is claimed before the archive query and positional validation finish.** Existing literature reports NIRCam multi-epoch observations from PRIMER, COSMOS-Web, COSMOS-3D; these are the known starting point, not an exhaustive inventory.

## Coordinates
`target.json` currently contains RA=150.0634341 deg, Dec=+2.3552216 deg (provisional coordinates carried over from our earlier discussion). **Verify against the authoritative COSMOS2025 entry or paper astrometry before final publication-level interpretation.** This matters for footprint boundaries. Search radius is used for discovery, not for asserting coverage.

## How to run
```bash
python -m pip install -r requirements.txt
python scripts/search_mast_epochs.py --target target.json --out results
```
Requires internet access to MAST. The GitHub Actions workflow runs on push / manually / hourly, and commits the *metadata-only* output under `results/`. No proprietary FITS files, images or large data are downloaded.

## Outputs
- `results/REPORT.md`: readable summary, status and warnings.
- `results/observations.csv`: individual MAST observation records with footprint status.
- `results/epoch_groups.csv`: grouped dates, missions, instruments and filters. **Candidate epoch**, not necessarily a verified astrometric epoch.
- `results/science_products.csv`: relevant science FITS file metadata for follow-up validation.
- `results/status.json`: machine-readable execution state and source.
- `results/other_missions.csv`: nearby records from non-HST/JWST archives, as context.

## Scientific distinction: archive overlap is NOT a confirmed detection
1. MAST positional discovery provides candidate observations. Spatial footprint (`s_region`) intersection is checked when possible.
2. Observation footprint can include detector gaps; a **pixel-level verification with the calibrated FITS science footprint/gWCS** is still needed.
3. A covered exposure is not yet a **detection of Shabgard**. Require PSF/centroid and SNR assessment in that filter.
4. Individual exposures/dithers on one visit/night **do not** create independent parallax phases.
5. HST optical non-detection is useful for color limits, but is **not** an astrometric epoch unless the target is detected.
6. Cross-instrument JWST NIRCam/MIRI or HST/JWST astrometry must characterize PSF, distortion, chromatic registration, and independent galaxy controls.

## Next-stage scientific audit
For each candidate exposure in the results, check `s_region` against its specific footprint, then download only necessary headers/small cutouts to test gWCS/pixel coverage and source S/N. Extract DATE-OBS, TIME-OBS, MJD, filter, exact RA/Dec, and local stationary-reference quality. Build a per-epoch table with status `covered`, `detected`, `centroidable`, and predicted parallax factor using JWST or HST ephemerides. Only then evaluate five-parameter PM + parallax identifiability and injection-recovery.

Original publication: https://arxiv.org/abs/2610.07861
MAST archive: https://mast.stsci.edu/
