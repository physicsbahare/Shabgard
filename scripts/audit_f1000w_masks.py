"""Diagnose rejected F1000W apertures without relaxing quality masks."""
import glob,json
import numpy as np
from astropy.io import fits
from jwst import datamodels
from miri_calibration_pilot import sky_samples, exposure_mjd
from pathlib import Path
rows=[]
for name in glob.glob("miri_data/*_mirimage_cal.fits"):
    with fits.open(name) as f, datamodels.open(name) as m:
        if f[0].header.get("FILTER")!="F1000W": continue
        a,e,q=[f[k].data for k in ("SCI","ERR","DQ")]
        mjd=exposure_mjd(f[0].header,f["SCI"].header)
        yy,xx=np.indices(a.shape)
        for t,x,y in sky_samples(mjd,m.meta.wcs):
            if not np.isfinite(x+y): continue
            r=np.hypot(xx-x,yy-y)
            for region,mask in [("core",r<=2),("aperture",r<=3.5),("annulus",(r>=7.35)&(r<=14.5))]:
                flags,counts=np.unique(q[mask],return_counts=True)
                good=np.isfinite(a)&np.isfinite(e)&(e>0)&((q&1)==0)
                rows.append(dict(file=Path(name).name,mjd=mjd,ref=t,x=x,y=y,region=region,n=int(mask.sum()),good=int((mask&good).sum()),nonfinite_sci=int((mask&~np.isfinite(a)).sum()),nonfinite_err=int((mask&~np.isfinite(e)).sum()),nonpositive_err=int((mask&(e<=0)).sum()),dq={str(int(k)):int(v) for k,v in zip(flags,counts)}))
Path("miri_results/f1000w_mask_audit.json").write_text(json.dumps(rows,indent=2))
print(json.dumps(rows,indent=2))
assert rows,"No F1000W detector samples"
