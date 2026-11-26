# Free-centroid MIRI photometry diagnostic summary, 2026-11-26
Recomputed directly from the 21 per-exposure JSON records printed by successful run 45553291617, job 138999273864. All values in uJy except recovery fraction and PSF sum. Sample SD across dithers, not standard error; null noise from fixed-position fit is not calibrated free-search uncertainty. Stress-test median flux divides by median recovery from artificial broadened PSF only, NOT recommended correction or upper bound. PSF sum is finite template normalization, not a measured aperture correction. Masked F1000W exposures absent, not nondetections. No final flux or SED fitting.

```json
[
  {
    "filter": "F770W",
    "n": 8,
    "median": 0.9089821925201587,
    "sd": 0.06041635429613942,
    "mad": 0.05856108479019033,
    "range": [
      0.8531301475176462,
      1.0364686220461972
    ],
    "formal_median": 0.05893854645275903,
    "null_proxy_median": 0.08543400921391392,
    "recovery": 0.8393852806886422,
    "stress_test_median_flux": 1.0829141437581777,
    "psf_sum": 0.9488309311427038,
    "controls": 96,
    "exceedances": 0
  },
  {
    "filter": "F1000W",
    "n": 1,
    "median": 2.116671988613333,
    "sd": null,
    "mad": 0,
    "range": [
      2.116671988613333,
      2.116671988613333
    ],
    "formal_median": 0.13498319517965524,
    "null_proxy_median": 0.18743867319965085,
    "recovery": 0.8402916784981924,
    "stress_test_median_flux": 2.518972926634649,
    "psf_sum": 0.9349064374612487,
    "controls": 12,
    "exceedances": 0
  },
  {
    "filter": "F1280W",
    "n": 6,
    "median": 2.241868454291371,
    "sd": 0.2813581190350731,
    "mad": 0.17869269979136057,
    "range": [
      1.7776742741822313,
      2.4579125058379754
    ],
    "formal_median": 0.22830184125966177,
    "null_proxy_median": 0.33524180308106466,
    "recovery": 0.8516500287602866,
    "stress_test_median_flux": 2.6323822915320867,
    "psf_sum": 0.9298623838505062,
    "controls": 72,
    "exceedances": 0
  },
  {
    "filter": "F1500W",
    "n": 6,
    "median": 2.2029883472653293,
    "sd": 0.27418230558569845,
    "mad": 0.24973986926724076,
    "range": [
      1.9871255452323056,
      2.678220249554524
    ],
    "formal_median": 0.3207542165503993,
    "null_proxy_median": 0.3240326436897088,
    "recovery": 0.8085794016966272,
    "stress_test_median_flux": 2.7245170265812355,
    "psf_sum": 0.921556608878153,
    "controls": 72,
    "exceedances": 0
  }
]
```

Next: retain full matched-search null and injection samples, vary PSF width and sky, use empirical unresolved calibrators and joint covariance; compare F770W catalog 0.95145434 uJy with image association. Target relative astrometry run 45745600504 initiated using original CAL gWCS, matched centroid and registration artifact; no absolute tie yet.
