# Independent fifth epoch execution — 2026-11-26

Actual successful fourth-epoch job 139869205481 (run 45808293268): parsed FOURTH_EPOCH and FACTORS JSON. 22 common baseline candidate references; 11,11,13 retained at scatter cut 10 mas. Three distinct F140M CAL centroid results east -181.7823,-173.0965,-166.2625 mas; north -124.1744,-125.9789,-119.6452 mas. Sample SD [7.778321052715062,3.2630375031536243] mas; not final epoch covariance. No stationary classification verified.

Geometry-only screen: actual four baseline median MJD and rounded JWST Horizons factors [[60082.30187,-0.96338,0.34522],[60308.94067,0.70215,-0.36309],[60783.83768,-0.83896,0.21836],[61183.22692,-0.962786,0.36276]]. For each hypothetical fifth date, use the same winter factor (+0.70215,-0.36309) as a proxy, not candidate ephemerides. Weighted sum of squared parallax-factor residuals after projecting out intercept and time in each coordinate; baseline 2.0219528736643593. Sigma ratios relative to baseline, assuming independent isotropic equal errors for baseline, and fifth errors equal, 3x, or 10x larger. No measured parallax or forecast of actual precision.
```json
[
  {
    "date": "2022-12-26",
    "ideal_sigma_ratio": 0.907909311611745,
    "sigma_ratio_if_3x_noisier": 0.9807151433703265,
    "sigma_ratio_if_10x_noisier": 0.9980693523190941
  },
  {
    "date": "2024-01-04",
    "ideal_sigma_ratio": 0.8245786661395059,
    "sigma_ratio_if_3x_noisier": 0.9675593376032169,
    "sigma_ratio_if_10x_noisier": 0.9968372068954191
  },
  {
    "date": "2024-12-21",
    "ideal_sigma_ratio": 0.7589784580878526,
    "sigma_ratio_if_3x_noisier": 0.9534623683203728,
    "sigma_ratio_if_10x_noisier": 0.995414649069324
  },
  {
    "date": "2025-12-29",
    "ideal_sigma_ratio": 0.7363948877537273,
    "sigma_ratio_if_3x_noisier": 0.9385115103581212,
    "sigma_ratio_if_10x_noisier": 0.9936482430443259
  },
  {
    "date": "2011-12-10",
    "ideal_sigma_ratio": 0.9760792220487615,
    "sigma_ratio_if_3x_noisier": 0.980789979246935,
    "sigma_ratio_if_10x_noisier": 0.9940707957180752
  },
  {
    "date": "2012-01-26",
    "ideal_sigma_ratio": 0.9764434455321691,
    "sigma_ratio_if_3x_noisier": 0.9811571274661118,
    "sigma_ratio_if_10x_noisier": 0.9942481184635605
  }
]
```
December 2025 highest geometric priority; December 2024 only one usable CAL; Jan 2024 close to already measured Dec 2023; HST longer proper-motion baseline but detection unverified. Actual old PRIMER F770W log 138502536750 reports 86 public candidate CALs and only 3 selected; prior sample does not rule out target coverage. New code eebc002cf4dedcf5cef05e13108858117455e52a maps previous physical STPSF target centroids to original CAL gWCS in F1280W/F1500W separately and measures direct reference centroids against the 22 NIRCam common controls, no assumption of achromatic morphology. Code 9fbc418ab0202f6fe3ec59803719ca2029ebe39d screens all 2022 F770W headers and actual pixels for plausible detector members (unknown-header downloads capped at 12). Workflow c09047248857d44ed84294a57de5e8a59e214617, run 45815669496 verified queued at launch; no result asserted. Next: inspect both actual job logs and numerical artifacts, compare F1280W/F1500W translations, residuals, common-reference leave-one-out and real Horizons factors; accept neither as validated without chromatic covariance. If 2022 pixels positive, fit physical PSF and assess same-filter tie to 2024 F770W; if insufficient, inspect HST WFC3/IR FLTs from 2010–2012 and compare predicted source with controls. Do not fit SED or claim parallax. No direct server access needed for these jobs.
