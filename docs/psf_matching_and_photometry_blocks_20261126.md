# Paired common-grid PSF matching and independent MIRI block repeatability — 2026-11-26

Astrometry: previous run 45760577380, job 139709332505 log confirmed 17 F140M measurements in 3 CALs and 24 F770W in 8 CALs, at 0.5 arcsec SNR proxy >5. Prior artifact 15089507016; previous numerical analysis in cross_instrument_aperture_sensitivity_20261126.md. Code commit 726fba9e3b399b5bd71d26ee82bb369adef29c2c runs paired original CAL common-grid centroids (0.04 arcsec), common intersection mask across all scenarios, same background annulus, raw vs Gaussian FWHM assumptions 0.23/0.26/0.30 arcsec and additional smoothing to 0.40 arcsec. F140M assumed 0.055 arcsec. Pair combinations share exposures; no independent-pair error bars. Equal galaxy medians, leave-one-galaxy-out and exclusion of 638656 only as sensitivity test. Gaussian approximation not a physical PSF or validated astrometric tie. Run 45774659681 confirmed queued when launched; no new numerical result asserted at checkpoint creation. Next: inspect this run's actual logs and psf_matched_summary.json plus stamps; compare exact same pairs and aperture spread; if successful implement physical PSF and morphological vetting, otherwise repair mask/grid issues.

Photometry: independently parsed all 21 actual per-exposure JSON output records from successful job 138999273864 (run 45553291617) and grouped by filter and filename exposure-sequence code. Blocks are descriptive, not independent epochs or calibrated errors. Results (uJy):
```json
{
  "stats": [
    {
      "block": "F770W_02201",
      "n": 4,
      "median_ujy": 0.955589189176814,
      "mean_ujy": 0.9680081763676475,
      "sd_ujy": 0.048387155674670626,
      "median_free_fixed_ratio": 2.0504253523132263,
      "max_null_target_ratio": 0.475601108409838
    },
    {
      "block": "F770W_04201",
      "n": 4,
      "median_ujy": 0.8748907170421552,
      "mean_ujy": 0.8741225653928774,
      "sd_ujy": 0.017301104646720192,
      "median_free_fixed_ratio": 1.759089007916799,
      "max_null_target_ratio": 0.6819376779275846
    },
    {
      "block": "F1000W_05201",
      "n": 1,
      "median_ujy": 2.116671988613333,
      "mean_ujy": 2.116671988613333,
      "sd_ujy": null,
      "median_free_fixed_ratio": 1.455130540610139,
      "max_null_target_ratio": 0.3388069493141677
    },
    {
      "block": "F1280W_02101",
      "n": 3,
      "median_ujy": 2.241995532900962,
      "mean_ujy": 2.250204820305968,
      "sd_ujy": 0.014439568690269555,
      "median_free_fixed_ratio": 1.0607947667413495,
      "max_null_target_ratio": 0.5035860847217444
    },
    {
      "block": "F1280W_04101",
      "n": 3,
      "median_ujy": 1.7830032387891086,
      "mean_ujy": 2.006196672936438,
      "sd_ujy": 0.3912064604993024,
      "median_free_fixed_ratio": 1.1005456827640248,
      "max_null_target_ratio": 0.5174403824507736
    },
    {
      "block": "F1500W_06101",
      "n": 3,
      "median_ujy": 2.0819566811110595,
      "mean_ujy": 2.208841155348995,
      "sd_ujy": 0.30559736614056776,
      "median_free_fixed_ratio": 1.1924931395037563,
      "max_null_target_ratio": 0.7559473381527058
    },
    {
      "block": "F1500W_08101",
      "n": 3,
      "median_ujy": 2.2455663563299115,
      "mean_ujy": 2.3613989813617278,
      "sd_ujy": 0.27765927652497513,
      "median_free_fixed_ratio": 1.1094469434588414,
      "max_null_target_ratio": 0.7399464381549815
    }
  ],
  "contrasts": [
    {
      "filter": "F770W",
      "block_median_difference_ujy": 0.08069847213465886,
      "fraction_of_grand_median": 0.08877893626378076
    },
    {
      "filter": "F1280W",
      "block_median_difference_ujy": 0.4589922941118536,
      "fraction_of_grand_median": 0.204736497020266
    },
    {
      "filter": "F1500W",
      "block_median_difference_ujy": -0.16360967521885206,
      "fraction_of_grand_median": -0.07426715416899425
    }
  ]
}
```
F770W catalog comparison reference 0.95145434 uJy, source association and cross-calibration still pending. No correction, significance test, or SED fit. F1000W only one usable exposure; masked exposures not nondetections. Next executable photometry action: compare block-specific field point-source flux ratios and background / detector position with target block contrast; derive covariance or systematic term rather than averaging away dither scatter.
