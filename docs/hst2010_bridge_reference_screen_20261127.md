Source: artifact 15169064227, downloaded and analyzed with Python. Cross-visit bridge residual norms for accepted references 10,36,41,43,6,64,92: 21.589,10.483,84.446,9.272,0.056,8.285,228.043 mas respectively. Small visit difference does not imply small absolute bridge residual (ref 92: 4.419 mas).

Diagnostic only: select accepted references by bridge residual norm threshold; compute per-FLT original translation minus unweighted median east/north of available selected references. Mean shift and sample SD across four FLTs (mas):

| threshold | refs | mean E,N | SD E,N |
|---|---|---|---|
|10|43,6,64|4.608,4.038|3.480,4.816|
|15|36,43,6,64|-0.836,7.858|3.875,1.836|
|25|10,36,43,6,64|4.199,8.186|3.647,5.610|
|100|10,36,41,43,6,64|-1.041,6.095|3.967,4.460|
|300|all seven|-6.280,8.186|7.946,5.610|

Reproduce: load registration.json[0], measurements.json, bridge_refs.json from ZIP; for each threshold select registration refs with numpy.linalg.norm(bridge_refs[id]['residual']) <= threshold; for each registration exposure select measurement rows matching filename and selected ref, take numpy.median of east,north axis=0, subtract from original translation; report mean and std(ddof=1) over exposures. Not a refit or independent validation; original estimator differs even for all seven. Do not use these shifts as calibrated covariance or target positions. Absolute residual screening is useful for flagging 41 and 92, not proof of stationary galaxy classification. Direct MAST retrieval of ibhm55e7q_flt.fits attempted in local runtime; DNS resolution failed. Pixel injection remains blocked locally; next use archive-capable workflow. No change to pilot parallax.