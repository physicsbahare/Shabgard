# Parallax validation plan
Baseline: six-epoch pilot 6.60 +/- 3.49 mas; nominal inverse 152 pc, not validated distance. See joint_six_epoch_pilot_20261127.md. Do not reduce errors by assumption.

1. Run physical/ePSF injection recovery on actual detector exposures, matching target flux, subpixel phase, background and masks; recover bias, 2D covariance and coverage, including shared registration uncertainty. Test stationary controls and false parallax.
2. Classify stationary compact galaxy references, test chromatic morphology, reject moving stars using evidence, leave-one-reference-out and cross-filter/instrument registration covariance.
3. Validate HST 2010 four FLTs, association, dither agreement and registration; include only if independent and calibrated. Use exact exposure midtimes and observatory barycentric vectors.
4. Search opposite-phase independent NIRCam epochs; quantify marginal Fisher leverage with full covariance. COSMOS-Web removal raises pilot sigma to 15.1 mas. Avoid redundant census and same-phase epochs.
5. Fit GLS five parameters versus PM-only, free signed parallax, full covariance, rank, residuals, leave-one-epoch-out and systematic sensitivity. Sample joint parameter distribution if warranted; no positivity or assumed distance prior for detection. Target sigma <=2.2 mas (3 sigma) or <=1.32 mas (5 sigma) conditional on pi=6.6, not acceptance by fiat.
6. Produce Liu et al Fig 3/4 style sky trajectory and PM-subtracted parallax ellipse using actual positions, exact factors and joint uncertainty samples; show epoch covariance ellipses, no fabricated exposure points or injection clouds. Label pilot until validated. Preserve numerical tables and reproducible scripts.
Maintain 70% astrometry/30% calibrated MIRI photometry; no SED modeling. Each hourly run must execute real work, verify run/artifact/commit, repair failures, and independently work if a job is pending. If data insufficient, document limits and feasible observations.
