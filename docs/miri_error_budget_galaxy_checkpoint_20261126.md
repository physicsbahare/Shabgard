# MIRI error-budget and image-based galaxy checkpoint, 2026-11-26

Actual ZIP artifact 15083964171 from successful run 45748539048 opened in local Python container. 1008 injection rows, 21 target fit rows; 252 unique exposure/angle controls. Only 364 injection rows pass both optimizer-success and nonboundary criteria for null and injected fits (F770W 156, F1000W 12, F1280W 128, F1500W 68). This is 91 distinct control positions after flags; same positions repeated at four widths, not independent. Strict subset selection itself can bias results. No target flux correction endorsed.

| Filter | N exp | Median uJy | Dither SD uJy | Median formal sigma uJy | Matched free-search null MAD uJy | Strict controls/total |
|---|---:|---:|---:|---:|---:|---:|
| F770W |8|0.908982|0.060416|0.058939|0.068648|39/96|
| F1000W |1|2.116672|undefined|0.134983|0.132742|3/12|
| F1280W |6|2.241869|0.281358|0.228302|0.294380|32/72|
| F1500W |6|2.202990|0.274183|0.320754|0.398892|17/72|

Median recovery at blur widths 0,0.25,0.45,0.65 pix (all vs strict):
F770W all 0.901992,0.901677,0.835323,0.724011; strict 0.911814,0.911499,0.845336,0.735570.
F1000W all 0.886743,0.886539,0.842993,0.766130; strict 0.995243,0.995064,0.956999,0.885541.
F1280W all 0.877982,0.877847,0.849331,0.798315; strict 0.925289,0.925162,0.897662,0.846851.
F1500W all 0.821419,0.821311,0.798702,0.760728; strict 0.896321,0.896225,0.876207,0.839539.
Null MAD includes positive search bias and local background/confusion; not an unbiased Gaussian target error. Zero control exceedances do not establish calibrated FAP. Reproducible script scripts/miri_error_budget.py commit db84886587e87043590f3e357c56810befa2211d, derived numerical statistics above executed independently in local Python. F1000W masked exposures excluded. No final flux uncertainties or SED fitting.

Next executable stage: direct image centroids on original F770W eight CAL and F140M three CAL for six candidate galaxy IDs, matched angular apertures 0.35/0.5/0.7 arcsec, local robust plane and signed moments, diagnostic stamps. Code b7f94e68bdc86bd7fce7c5268cb05740ad020055 and workflow 94650f09c18a90ebc80fce5795f90e658b2fc989. Not yet verified execution result. Morphology, PSF/chromatic offsets, nearby blends, absolute registration and target covariance still require validation. Do not treat apertures as independent or use a catalog-only offset for parallax.
