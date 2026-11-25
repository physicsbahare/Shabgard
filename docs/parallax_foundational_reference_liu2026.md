# Foundational reference for the Shabgard parallax project

Identified by the project lead on 2026-11-25 as the paper that inspired this project.

## Reference

Liu, Feng-Yuan Frey, et al. (2026), *The cold substellar frontier in JWST Extragalactic surveys: a robust parallax measurement yields the most distant known ultra-cold Y dwarf*. Submitted to MNRAS; arXiv:2611.08346v1, 5 November 2026.

- [HTML, version 1](https://arxiv.org/html/2611.08346v1)
- [Abstract](https://arxiv.org/abs/2611.08346v1)
- [PDF](https://arxiv.org/pdf/2611.08346v1)
- [DOI](https://doi.org/10.48550/arXiv.2611.08346)

## Why this is the benchmark

The preprint reports a geometric distance to JADES-GS-BD-5 from calibrated multi-epoch NIRCam astrometry over approximately two years: parallax 15.7 (+1.4/-1.3) mas, distance 64 (+6/-5) pc, and total proper motion 49.8 +/- 0.3 mas/yr. These are the published preprint's values, not Shabgard measurements.

## Shabgard implementation plan

Target: COSMOS2025 ID 639188. Use the Liu et al. preprint as the primary methodological reference for the astrometric track. The following are project requirements, not a claim that every item reproduces the paper exactly:

1. Recover original exposure-level positions, midtimes, uncertainties and stationary reference sources from the existing Shabgard proper-motion analysis; verify the effective COSMOS2025 B7 coordinate epoch rather than fixing the provisional MJD 60308.
2. Validate the new 2026-05-18 NIRCam F140M centroid measurements using physical/empirical PSFs and stationary-field registration. Include MIRI visits only after source association, centroid validation and a quantified cross-instrument reference tie. Include HST only if a usable detection and registration are established.
3. Establish a common astrometric reference frame, per-exposure transformations, realistic centroid and registration covariance, chromatic and instrument-dependent offsets; test leave-one-epoch-out and cross-filter consistency.
4. Compute observer-specific parallax factors at the actual exposure midtimes using JWST/HST ephemerides. Check cadence and identifiability before a five-parameter fit of reference position, two proper-motion components and parallax. Compare with a proper-motion-only solution and validate the false-parallax rate using controls/injections.
5. Document measured versus candidate epochs, residuals, uncertainties and reproducible code. Do not claim a geometric distance until the parallax is validated.

The atmospheric analysis is separate: complete calibrated MIRI photometry and its uncertainty budget, but stop before model SED fitting until the project lead provides the server model grids.

## Pending inputs from project lead

- Original astrometry epoch table / analysis code / reference-source files.
- COSMOS-Web B7 SCI/WHT mosaic and exposure-contributor information or server paths.

## Verification status

Bibliographic details and reported results checked against the preprint abstract via an accessible index. Direct HTML retrieval was unavailable during this update; the detailed method and appendices still need to be read and mapped against the implementation before claiming methodological reproduction. No new astrometry result is asserted here.

Related: [MIRI calibrated-detector checkpoint](miri_cal_science_checkpoint_20261124.md).
