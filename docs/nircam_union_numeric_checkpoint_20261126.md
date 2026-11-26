# Expanded F115W reference registration
Run 45795708087, job 139825682629, artifact 15105377048. Parsed actual REGISTRATION JSON from job log and computed per-program medians and sample SD in tangent-plane mas relative to catalog coordinate (not a catalog epoch measurement). All 11 exposures now have >=8 references; stationary classification remains false. First exposure is reference anchor, its near-zero scatter is not independent precision. No parallax fit.
```json
[
  {
    "program": "1837",
    "n": 4,
    "refs": [
      19,
      20,
      9,
      11
    ],
    "east": 35.10340210424653,
    "north": 62.29409129412389,
    "sd_e": 2.0109116219400485,
    "sd_n": 1.1303756224985169,
    "mjd": 60082.30187090856
  },
  {
    "program": "1727",
    "n": 4,
    "refs": [
      12,
      8,
      13,
      18
    ],
    "east": 2.427977200123331,
    "north": 17.703555448012054,
    "sd_e": 1.1386070029553577,
    "sd_n": 0.9981779726266644,
    "mjd": 60308.94067113426
  },
  {
    "program": "5893",
    "n": 3,
    "refs": [
      14,
      9,
      17
    ],
    "east": -97.33740721856817,
    "north": -57.59262121287634,
    "sd_e": 2.019787706197193,
    "sd_n": 1.995729785027538,
    "mjd": 60783.83768326389
  }
]
```
Next executable action: use reference measurement CSV from artifact to solve joint source positions and exposure translations, leave-one-source-out residuals and galaxy classification; validate target PSF and observer factors before fitting. Do not treat within-visit scatter as total uncertainty.
