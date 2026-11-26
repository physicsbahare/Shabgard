"""F770W cross-visit target displacement using actual CAL gWCS and network covariance."""
import csv,json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
from astroquery.jplhorizons import Horizons
from jwst import datamodels
from audit_calibrated_detector_local import exposure_mjd,RA,DEC
OUT=Path('displacement');OUT.mkdir(exist_ok=True)
def file(root,name):
    q=list(Path(root).rglob(name))
    if len(q)!=1:raise ValueError((name,q))
    return q[0]
def read(root,name):
    with file(root,name).open() as f:return list(csv.DictReader(f))
def main():
    center=SkyCoord(RA,DEC,unit='deg')
    diag={q['filename']:q for q in read('primer','pixel_audit.csv')}
    psf={q['filename']:q for q in read('network','miri_stpsf_diagnostic.csv')}
    old={q['filename']:q for q in read('positions2024','target_relative_positions.csv')}
    net=json.loads(file('network','network.json').read_text())
    rows=[]
    for n,d in diag.items():
        if n not in psf:continue
        q=psf[n]
        if q['free_optimizer_success'].lower()!='true' or q['free_near_search_boundary'].lower()=='true':raise ValueError(n)
        x=float(d['x_photometry'])+float(q['free_dx_pix']);y=float(d['y_photometry'])+float(q['free_dy_pix'])
        with fits.open(Path('cal')/n,memmap=False) as h,datamodels.open(str(Path('cal')/n)) as m:
            ra,dec=map(float,m.meta.wcs(x,y));pos=SkyCoord(ra,dec,unit='deg')
            a,b=center.spherical_offsets_to(pos);t=exposure_mjd(h[0].header,h['SCI'].header)
            rows.append(dict(filename=n,visit='2022',mid_mjd=t,east_mas=a.mas,north_mas=b.mas,flux_ujy=float(q['free_flux_ujy'])))
    for n,q in old.items():
        rows.append(dict(filename=n,visit='2024',mid_mjd=float(q['mid_mjd']),east_mas=float(q['raw_east_mas']),north_mas=float(q['raw_north_mas']),flux_ujy=float(q['free_flux_ujy'])))
    assert len(rows)==12 and len(set(r['filename'] for r in rows))==12
    lookup={q['filename']:q for q in rows}
    results=[]
    for k in net:
        if 'error' in k:continue
        off={q['filename']:q for q in k['offsets']}
        if set(off)!=set(lookup):continue
        groups={}
        for n,q in lookup.items():
            v=np.array([q['east_mas'],q['north_mas']])-np.array([off[n]['ra_mas'],off[n]['dec_mas']])
            groups.setdefault(q['visit'],[]).append(v)
        a,b=np.array(groups['2022']),np.array(groups['2024'])
        res=dict(radius_mas=k['radius_mas'],snr_cut=k['snr_cut'],displacement_2024_minus_2022_mas=(b.mean(axis=0)-a.mean(axis=0)).tolist(),median_displacement_mas=(np.median(b,axis=0)-np.median(a,axis=0)).tolist(),scatter_2022_mas=a.std(axis=0,ddof=1).tolist(),scatter_2024_mas=b.std(axis=0,ddof=1).tolist())
        results.append(res)
    cov=np.load(file('network','relative_translation_covariance.npz'),allow_pickle=False)
    names=list(cov['filenames']);assert set(names)==set(lookup)
    weights=np.array([-1/4 if lookup[n]['visit']=='2022' else 1/8 for n in names])
    A=np.kron(weights.reshape(1,-1),np.eye(2));C=A@cov['covariance_mas2']@A.T
    times=Time(sorted(q['mid_mjd'] for q in rows),format='mjd',scale='utc')
    vec=Horizons(id='-170',location='500@0',epochs=times.tdb.jd.tolist()).vectors(refplane='earth')
    if not np.allclose(vec['datetime_jd'],times.tdb.jd,rtol=0,atol=1e-7):raise ValueError('ephemeris mismatch')
    xyz=np.array([vec[k].quantity.to_value(u.au) for k in ('x','y','z')]).T
    ra,dec=np.deg2rad([RA,DEC]);east=np.array([-np.sin(ra),np.cos(ra),0]);north=np.array([-np.cos(ra)*np.sin(dec),-np.sin(ra)*np.sin(dec),np.cos(dec)])
    factors={float(t):[-float(v@east),-float(v@north)] for t,v in zip(times.mjd,xyz)}
    for q in rows:q['px'],q['py']=factors[q['mid_mjd']]
    with (OUT/'positions.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary=dict(results=results,translation_only_contrast_covariance_mas2=C.tolist(),stationary_verified=False,parallax_ready=False)
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2));print('DISPLACEMENT',json.dumps(summary),flush=True)
    if not results:raise RuntimeError('No usable network solution')
if __name__=='__main__':main()
