# F770W all-dither relative registration checkpoint (2026-11-25)

Input: run [45582311047](https://github.com/physicsbahare/Shabgard/actions/runs/45582311047), artifact [15006733360](https://github.com/physicsbahare/Shabgard/actions/runs/45582311047/artifacts/15006733360), field_results/field_candidates.csv, 181 actual CAL-derived centroid rows. ZIP downloaded and opened in the current execution environment. Script: scripts/miri_f770_network_registration.py, commit 99371c96b75826364493d401a80428ebdd406a1d. Executed locally with NumPy/SciPy; initial 300-bootstrap version timed out, reran with 60 bootstrap draws and 5 IRLS iterations. 
Method: mutual nearest neighbors between every exposure pair in tangent plane, unique-source graph groups (reject groups containing duplicate exposure IDs), joint per-source positions and per-exposure translations with a fixed gauge and robust radial reweighting. Source-block bootstrap and leave-one-source-out tests; no claim of independent pairwise edges. All radii/SNR thresholds diagnostic, not optimized or calibrated.

| Match radius (mas) | SNR proxy cut | Unique groups | Rows | Pair links | Residual MAD RA*,Dec (mas) | Valid bootstrap / 60 | Valid LOO |
|---|---|---|---|---|---|---|---|
|500|0|74|169|154|10.47,12.96|59|74|
|250|0|73|167|153|10.74,13.08|58|73|
|150|0|72|165|152|10.63,12.63|57|72|
|250|15|12|45|92|5.86,4.93|58|12|
|250|30|4|20|50|5.69,8.77|53|4|

All 28 exposure pairs have at least one matched candidate for all trials. All eight exposure translations identifiable in full sample. Reference is jw01727016001_02201_00003_mirimage_cal.fits in all trials. At 250 mas / all SNR, offsets (measured minus common frame; subtract for relative registration):

| Exposure suffix | Distinct controls | RA* offset (mas) | Dec offset (mas) | Bootstrap SD RA*,Dec (mas) | Max LOO shift (mas) |
|---|---|---|---|---|---|
|02201_00001|13|-12.94|18.00|9.10,6.01|4.68|
|02201_00002|11|-4.24|9.53|9.49,7.14|5.61|
|02201_00003|33|0|0|0,0 (gauge)|0|
|02201_00004|26|-5.07|5.82|9.16,7.60|11.13|
|04201_00001|13|-4.21|16.40|8.75,5.31|4.84|
|04201_00002|12|-13.95|12.12|9.90,5.50|6.42|
|04201_00003|33|-0.19|-6.34|3.74,5.74|1.53|
|04201_00004|26|-6.71|10.98|10.88,9.38|11.35|

Only 4 unique groups remain at SNR>30 (threshold implemented >=30); the relative Dec offset of 04201_00003 shifts from -6.34 to -32.77 mas and bootstrap SD is 33.44 mas. Therefore small residuals in bright subsets are not proof of accuracy. Residual MAD is within-source post-fit scatter, not absolute or predictive astrometric error, and not directly comparable to prior single-pair residuals. Bootstrap conditional on connectivity and small number of draws; not full covariance or final error budget. No stationary classification, no cross-epoch tie, no target correction applied, no parallax measurement. 

Next executable milestone: use matched-group sky positions to identify stationary extragalactic sources in NIRCam/COSMOS2025; recover original CAL target centroid gWCS coordinates and apply relative corrections with full covariance (including common reference uncertainty and target centroid errors). Validate cross-instrument morphology/PSF and association; only then consider F770W as an astrometric epoch. Do not fit SED.
