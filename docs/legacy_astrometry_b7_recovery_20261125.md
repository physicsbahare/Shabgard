# Shabgard legacy astrometry and COSMOS-Web B7 recovery checkpoint — 2026-11-25

Read this before continuing the hourly parallax work, together with [coordinate epoch policy](shabgard_coordinate_epoch_working_policy_20261121.md), [MIRI checkpoint](miri_cal_science_checkpoint_20261124.md) and [foundational reference](parallax_foundational_reference_liu2026.md).

## Recovered original assets
- [Legacy notebook](../data/legacy_astrometry/proper_motion_fig4_corrected_baseline.ipynb): nine original code cells and saved textual outputs from the August 2026 notebook recovered from the conversation file library. Embedded PNG outputs removed for portability; not a fresh run. Original notebook SHA256 125495a92d9e415bd80d7bd7242f00e3a5b27b1ec0a79160f68a19e16de77ccc.
- [Original baseline patch](../data/legacy_astrometry/proper_motion_baseline_fig4_patch.py): exact original file (git blob 002d7d00010c49802d31286ebc1c1860713f5434). Requires existing notebook helpers, reads original i2d headers and uses a FIXED assumed distance 160 pc. Do not use its parallax-subtracted positions as independent measurements of new parallax.
- [Recovered pairwise CSV](../data/legacy_astrometry/shabgard_recovered_original_pairwise_pm.csv): transcription of the stored notebook table, five rows, 18 columns, rounded to displayed precision; contains raw and reference-offset corrected displacements, reference scatters and preliminary PM using old calendar-date baselines. Not the original image-level CSV, not independent per-epoch RA/Dec, not final paper PM, not independent likelihood rows (shared epochs). Reference scatters are not centroid covariance.
- [Server inventory](../data/legacy_astrometry/targeted_inventory.txt): exact read-only inventory from 2026-11-23 (git blob 4dd2651859a7a44fc72a3253f9d1230f7241ab20). It is a limited inventory, not a complete exposure census.

Target COSMOS2025 639188 / COSMOS-3D 17985; catalog RA=150.06343413982992, Dec=+2.355221569555261 deg. Catalog epoch is unverified. Approx PM (-70,-54) mas/yr for searching only. Prior conversation reports published total PM 87.3 +/- 0.4 mas/yr; verify against paper before use. Earlier notebook pairs use 1.968515 yr; historical conversation reports corrected MJD approximately 60078.737 and 60782.536 and corrected baseline ~1.9269 yr; these rounded numbers are leads only, not a validated per-filter epoch table. Recover original header audit.

## Association-level i2d products in legacy code
- PRIMER PID 1837: jw01837-o003_t003_nircam_clear-f115w_i2d.fits, jw01837-o003_t003_nircam_clear-f200w_i2d.fits, jw01837-o003_t003_nircam_clear-f356w_i2d.fits
- COSMOS-Web PID 1727: jw01727-o051_t051_nircam_clear-f115w_i2d.fits
- COSMOS-3D PID 5893: jw05893-o014_t014_nircam_clear-f115w_i2d.fits, jw05893-o014_t014_nircam_clear-f200w_i2d.fits, jw05893-o014_t014_nircam_clear-f356w_i2d.fits

Retrieve MJD-AVG/EXPMID (or explicit start/end fallback with provenance), JWST observer position keywords and constituent exposure metadata from original primary headers/association products. Old cutouts dropped non-WCS metadata. Association-level dates do not by themselves establish detector-level coverage or catalog effective epoch.

## User workstation recovery
Historical notebook working directory lead: /home/bahareh/Desktop/Important_Codes/ (not verified live). Search subfolders proper_motion_image_level/, proper_motion_corrected_baseline/, ID639188_fake_source_PM_test/, jwst_image_cutouts_for_pm/, ID639188_Fig3_style_PM/.
Priority files from original code:
- image_level_pm_summary.csv
- ID639188_original_i2d_header_epoch_audit.csv
- ID639188_corrected_pair_values.csv
- ID639188_fake_source_trials.csv
- ID639188_fake_source_recovery_summary.csv
- image_level_pm_summary_corrected_baseline_parallax.csv
- ID639188_fake_source_trials_corrected_baseline_parallax.csv
- ID639188_fake_source_recovery_summary_corrected_baseline_parallax.csv
- ID639188_fake_source_combined_pm_samples_corrected.csv
- ID639188_fake_source_combined_pm_summary_corrected.csv
Seek per-image target and reference-source centroids (sky and pixel), IDs, quality flags, covariance, reference-frame offsets and injection results. The original independent per-epoch coordinate/error table and individual galaxy reference positions are NOT yet recovered.

## B7 on NTHU server
Case-sensitive path confirmed by saved inventory: /mnt/E_202/JWST/COSMOS-Web/NIRCam_mosaics/
Both COSMOS-Web and COSMOS_Web exist; use the hyphenated path for these mosaics. Four confirmed F444W files (prefix mosaic_nircam_f444w_COSMOS-Web_30mas_B7_v1.0):
- _sci.fits.gz (1.6507 GB)
- _wht.fits.gz (1.4986 GB)
- _err.fits.gz (1.4431 GB)
- _i2d.fits.gz (9.0813 GB)
Catalog /mnt/E_202/JWST/COSMOS-Web/COSMOSWeb_mastercatalog_v1.1.fits. Auxiliary directory /mnt/E_202/JWST/COSMOS-Web/auxiliary. Other bands and association/context/exposure provenance not yet verified.
User ran find on this root with stdout redirected to ~/shabgard_B7_inventory.txt on 2026-11-25; no terminal output is normal, not evidence of zero results. Ask user to run wc -l ~/shabgard_B7_inventory.txt and grep -i 'B7' ~/shabgard_B7_inventory.txt | head -40, then share inventory. No direct server filesystem access in this chat.

## Next execution priorities
1. Do not stall on missing workstation files: retrieve original association i2d metadata via MAST and identify calibrated exposures, reconstruct target/control centroid data if necessary. Keep original and reconstructed data distinct; avoid repeating work once verified.
2. Inspect B7 association/CON/context and local exposure weights at the target for all four detection bands; determine effective date(s) and parallax factors of the catalog detection measurement, not release year, not F444W WHT alone. If local file unavailable use archive provenance and state limitation.
3. Register to stationary galaxies, retain per-epoch centroid errors/covariance, use actual observer position and mid-exposure times; include F140M and validated MIRI positions only after source association and cross-instrument systematics checks.
4. Fit position, PM and parallax jointly only with adequate independent epochs and parallax-factor diversity. Never fix 160 pc for the new fit or treat old pairwise rows as independent. Continue the MIRI photometric calibration track in parallel as specified in the MIRI checkpoint.
