# Six-epoch joint astrometry pilot — 2026-11-27
Executed locally with Python on downloaded actual artifacts: baseline NIRCam 15105377048, F140M 15110377673, HST refined 15150549729. Code commit 015d7ba3441f5458446887b4c88e7905562f6954. Positions below are [MJD, P_alpha, P_delta, east_mas, north_mas] in the pilot tangent plane:
```json
[[55905.822315288744,0.8667412415663365,-0.37214112491901014,850.3147373105685,761.1731371073772],
[55952.399801412495,0.34959847164242785,-0.28242362027230367,831.3205960475755,748.386631145115],
[60082.30187090856,-0.9633827831306112,0.3452214683902854,35.246668717064054,62.18449023309702],
[60308.94067113426,0.7021291437756914,-0.36309050214770267,2.571243813654974,17.59395438618583],
[60783.83768326389,-0.8389585778102971,0.2183614399171372,-97.1941406028562,-57.7022222739032],
[61183.22691935185,-0.9627863668938064,0.3627603946543518,-173.0965079697645,-124.17438772978367]]
```
HST factors are **approximate heliocentric Earth** from Swiss Ephemeris Moshier at individual midtimes, NOT actual HST barycentric factors. JWST factors are existing exposure-level Horizons values. HST positions reconstruct the exact pilot search prior plus common3 combined-filter peak offsets, equal filter weighting. 2010 HST excluded (unregistered); MIRI differential not included. Covariance uses HST dither mean covariance, independent hypothetical epoch floors, rank-one correlated HST common3-minus-common4 shifts [-2.5,-10,0,-15] mas and optional common HST frame covariance. Full rank 5. No positivity/distance prior. Results [HST floor, F115W floor, F140M floor, shared HST floor] mas:
- [5,3,8,0]: pi=6.598594 +/-2.105487 mas, delta chi2=9.821952, chi2=0.943004
- [10,5,10,10]: pi=6.600406 +/-3.490106 mas, delta chi2=3.576549, chi2=0.408391
- [20,5,15,20]: pi=6.578377 +/-3.528319 mas, delta chi2=3.476177, chi2=0.193502
- [30,8,20,40]: pi=6.577660 +/-5.633767 mas, delta chi2=1.363155, chi2=0.082864
Moderate PM east/north = -69.833784,-61.524986 mas/yr. Moderate leave-one-epoch-out [pi,sigma] mas: HST2011 [6.4950,3.5018], HST2012 [6.4810,3.5202], PRIMER [6.0493,4.3193], COSMOS-Web [11.8423,15.1192], COSMOS-3D [6.7415,3.6684], F140M [6.6955,3.5316]. Free common HST frame offset gives 6.703592 +/-3.514534 mas. Very small chi2 does not establish calibrated covariance; systematic floors are illustrative, not measured. COSMOS-Web provides critical phase leverage. No validated parallax or distance. Gaia audit run 45910179982 failed with TAP syntax error on all queries; repaired script commit e8d9af3d23e7abc6947ecf51cefc9dc7807ff446 and trigger commit a1b9f816d481b6bffd1d8ff326b22883a5532d47; rerun 45915718672 verified in progress. Next: obtain actual HST Horizons vectors, confirm source association/stationary references, NIRCam cross-filter systematics, estimate realistic covariance and investigate alternative winter epoch. No SED fitting.
