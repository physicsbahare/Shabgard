# Cross-instrument aperture registration sensitivity — 2026-11-26
Run 45760577380 job 139709332505 succeeded; inspected actual logs and downloaded ZIP 15089507016, read galaxy_results/galaxy_aperture_centroids.csv locally with pandas. At radius 0.5 arcsec and SNR proxy >5, 17 F140M and 24 F770W measurements across 3 and 8 exposures respectively. Equal weight per galaxy: compute mean east/north offset per filter then subtract F140M from F770W. No stationary classification or final frame tie.
| radius arcsec | all six median east,north mas | excluding 638656 median east,north mas | five-source sample SD east,north mas |
|---|---|---|---|
|0.35|49.34,50.07|63.87,56.38|43.11,20.61|
|0.50|56.09,61.80|82.14,63.43|47.97,17.35|
|0.70|103.92,76.90|132.49,86.83|34.05,45.43|
At 0.5 arcsec individual differences (east,north mas): 638656 (-83.82,-318.27), 638807 (91.91,63.43), 639322 (120.97,79.70), 639535 (30.04,60.18), 639727 (3.07,93.29), 639977 (82.14,49.23). Excluding 638656 is a sensitivity test, not an approved rejection. Five-source leave-one-out median ranges at 0.5 arcsec east 56.09–87.03, north 61.80–71.56 mas. Aperture dependence from 0.35 to 0.70 is 68.62 east,30.45 north mas after exclusion. This dominates a plausible few-mas parallax; do not use these translations as final astrometry. Unequal exposure coverage and chromatic morphology remain confounded. Next execute PSF matching of NIRCam to MIRI on same angular grid, use identical masks and compare per-galaxy offsets and leave-one-out; verify stationary morphology. No SED fitting.
