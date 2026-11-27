import csv,json,zipfile,sys
from pathlib import Path
import numpy as np
from astropy.time import Time
from astroquery.jplhorizons import Horizons
root=Path(sys.argv[1])
with zipfile.ZipFile(root/'hst_refined_fit.zip') as z:
    e=json.load(z.open('exposures.json'))
mjd=np.unique([q['mid_mjd'] for q in e])
t=Time(mjd,format='mjd',scale='utc').tdb.jd
ra,dec=np.deg2rad([150.0634341,2.3552216])
E=np.array([-np.sin(ra),np.cos(ra),0]);N=np.array([-np.cos(ra)*np.sin(dec),-np.sin(ra)*np.sin(dec),np.cos(dec)])
def vectors(id):
    v=Horizons(id=id,location='@0',epochs=t.tolist(),id_type=None).vectors(refplane='earth',aberrations='geometric')
    assert len(v)==len(t)
    assert np.max(abs(np.asarray(v['datetime_jd'])-t))<1e-6
    return np.column_stack([np.asarray(v[k],float) for k in ['x','y','z']])
earth=vectors('399')
hst=vectors('-48')
d=np.linalg.norm(hst-earth,axis=1)
assert np.all((d>3e-5)&(d<6e-5)),d
with (root/'hst_factors.csv').open('w') as f:
    w=csv.writer(f);w.writerow(['mid_mjd','px','py','earth_px','earth_py','hst_earth_au'])
    w.writerows(zip(mjd,-hst@E,-hst@N,-earth@E,-earth@N,d))
print('HST_FACTORS', (root/'hst_factors.csv').read_text())
