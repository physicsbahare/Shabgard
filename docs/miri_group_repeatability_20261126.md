# Dither-group repeatability and active science runs — 2026-11-26
Input actual artifact 14994941067 of run 45553291617; downloaded ZIP, opened miri_stpsf_diagnostic.csv, recomputed group means from 21 successful fit rows using filename exposure-group component. No final photometry.
| Filter | Group | N | Mean uJy | Sample SD uJy |
|---|---|---:|---:|---:|
| F770W | 02201 | 4 | 0.96800818 | 0.048387 |
| F770W | 04201 | 4 | 0.87412257 | 0.017301 |
| F1000W | 05201 | 1 | 2.11667199 | undefined |
| F1280W | 02101 | 3 | 2.25020482 | 0.014440 |
| F1280W | 04101 | 3 | 2.00619667 | 0.391206 |
| F1500W | 06101 | 3 | 2.20884116 | 0.305597 |
| F1500W | 08101 | 3 | 2.36139898 | 0.277659 |

Absolute group-mean difference divided by mean of the two group means: F770W 10.1932%, F1280W 11.4655%, F1500W 6.6761%. Only two groups each, not sufficient to estimate a population group covariance. The spread may reflect noise, sky or detector/PSF systematics, not established variability. Do not shrink uncertainty by sqrt(N) without covariance validation. F770W median free flux is 4.4639% below verified catalog value 0.95145434 uJy, without proving source association or independent calibration. Maximum per-exposure matched-null maximum / target flux over each filter: 0.681938, 0.338807, 0.517440, 0.755947 respectively; not an empirical false-positive probability. Masked F1000W files are not nondetections.
Target gWCS run 45745600504 failed due to missing astroquery during CAL selection. Fixed dependency commit c8aed51ae20c4fe6d2e3109ff5db2e6134bf0a12; replacement run 45748501114 observed in progress, not yet verified science. Extended PSF-control commit a0c1a250b5df431814634da8c34d419a0d6aebb8 and workflow commit 4b1739bbb421405b6b46e85da8eb24176f0c4556, run 45748539048 observed in progress. New per-trial output includes null and injected flux, subpixel phase, blur width 0/0.25/0.45/0.65 px, optimizer/boundary flags, centroid errors. Reused blank locations are not independent. Next inspect actual results/artifacts of both runs and use controls for uncertainty sensitivity, verify external astrometric tie and target association. No SED fitting.
