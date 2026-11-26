# December 2022 positive pixel diagnostics — 2026-11-26
Run 45815669496 actual jobs 139894647310 (success) and 139894647659 (failure) logs read. 2022: 86 CAL products screened; 82 outside header polygon, four actual CALs in pointing jw01837004006, all good aperture/sky and positive signal. Artifact 15115906201 available; values below calculated directly from four PRIMER_PIXEL JSON log records, not from ZIP inspection.
```json
{
  "n": 4,
  "median_mjd": 59950.15722550926,
  "median_raw_aperture_ujy": 0.5702023758281214,
  "sample_sd_ujy": 0.02130453764852448,
  "snr_proxy": [
    11.588542369032062,
    9.700706475170454,
    10.066581582056326,
    9.90644857790758
  ],
  "blank_exceedances": [
    0,
    0,
    0,
    0
  ]
}
```
Uncorrected aperture flux, not source identification or final photometry. Blank apertures not independent trials. 2025: 12 target positions and actual Horizons factors at 61042.7288–61042.7907; 12 control measurements total, zero exposures registered; RuntimeError No cross-instrument registration. Artifact 15114742704. Do not use as validated epoch. Prioritize same-filter 2022/2024 F770W instead of repeated MIRI/NIRCam PSF detours. New workflow commit 51873c3b0f128abbe29f6a03c221206e61349988, run 45818321064 confirmed queued. Uses 2022 cutout artifact, physical STPSF target fitting, downloads four exact CALs for field controls and combines actual 2024 field-control artifact 15006733360 for same-filter network. Next inspect logs and numerical outputs, verify connectivity, cross-visit reference residuals and source association; if successful derive relative target displacement and uncertainty, tie to NIRCam only if warranted. No measured fifth epoch or parallax yet; no SED fit.
