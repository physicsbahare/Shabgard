"""Add May 2026 F140M to existing F115W frame; diagnostic, not validated parallax."""
import csv,json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
from jwst import datamodels
from nircam_same_filter_astrometry import download,fit_one,usable,offset
from audit_calibrated_detector_local import exposure_mjd,sky_at,RA,DEC
OUT=Path('fourth_epoch');OUT.mkdir(exist_ok=True)
def read(name):
    p=list(Path('baseline').rglob(name))
    if len(p)!=1:raise RuntimeError(f'Expected one {name}: {p}')
    return list(csv.DictReader(p[0].open()))
def write(name,rows):
    if rows:
        with (OUT/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=sorted(set().union(*(q.keys() for q in rows))));w.writeheader();w.writerows(rows)
def main():
    anchors={q['id']:q for q in read('anchor_refs.csv')}
    regs={q['filename']:q for q in read('registration_pilot.csv')}
    refs=read('reference_centroids.csv')
    by={}
    for q in refs:
        if q['filename'] not in regs:continue
        reg=regs[q['filename']]
        by.setdefault(q['ref_id'],{}).setdefault(q['program'],[]).append([float(q['east_mas'])-float(reg['east_mas']),float(q['north_mas'])-float(reg['north_mas'])])
    base={}
    for k,v in by.items():
        if len(v)!=3:continue
        a=np.array([np.median(x,axis=0) for x in v.values()])
        center=np.median(a,axis=0)
        p=SkyCoord(float(anchors[k]['ra']),float(anchors[k]['dec']),unit='deg').spherical_offsets_by(center[0]*u.mas,center[1]*u.mas)
        base[k]=(p,float(np.max(np.linalg.norm(a-center,axis=1))))
    print('BASELINE_COMMON',len(base),flush=True)
    matches=[];targets=[];summ=[]
    for n in ('00001','00002','00003'):
        name=f'jw07814003013_02201_{n}_nrca1_cal.fits'
        fn=download('mast:JWST/product/'+name,name)
        with fits.open(fn,memmap=False) as h,datamodels.open(fn) as m:
            a,e,d=h['SCI'].data,h['ERR'].data,h['DQ'].data;w=m.meta.wcs
            t=exposure_mjd(h[0].header,h['SCI'].header)
            q,p,_,_=fit_one(a,e,d,w,*sky_at(t,60308))
            target=dict(filename=name,mid_mjd=t,ra=p.ra.deg,dec=p.dec.deg,usable=bool(usable(q)),**q)
            targets.append(target)
            for k,(b,scatter) in base.items():
                try:
                    z,ps,_,_=fit_one(a,e,d,w,b.ra.deg,b.dec.deg)
                    if not usable(z,50):continue
                    dx,dy=offset(b,ps)
                    matches.append(dict(filename=name,ref_id=k,east=dx,north=dy,baseline_epoch_scatter_mas=scatter,**z))
                except Exception:continue
            rr=[x for x in matches if x['filename']==name]
            for threshold in (10,20,50,1e9):
                chosen=[x for x in rr if x['baseline_epoch_scatter_mas']<=threshold]
                if len(chosen)<3:continue
                v=np.array([[x['east'],x['north']] for x in chosen])
                c=np.median(v,axis=0);mad=np.maximum(1,1.4826*np.median(abs(v-c),axis=0))
                keep=np.all(abs(v-c)<=4*mad,axis=1);v=v[keep]
                if len(v)<3:continue
                c=np.median(v,axis=0);loo=np.array([np.median(np.delete(v,i,axis=0),axis=0) for i in range(len(v))])
                ps=p.spherical_offsets_by(-c[0]*u.mas,-c[1]*u.mas)
                dx,dy=offset(SkyCoord(RA,DEC,unit='deg'),ps)
                rec=dict(filename=name,mid_mjd=t,baseline_scatter_cut=threshold,n_ref=len(v),n_initial=len(chosen),target_usable=bool(usable(q)),registered_east_mas=dx,registered_north_mas=dy,translation_east_mas=float(c[0]),translation_north_mas=float(c[1]),reference_mad_mas=(1.4826*np.median(abs(v-c),axis=0)).tolist(),loo_range_mas=np.ptp(loo,axis=0).tolist(),stationary_verified=False)
                summ.append(rec);print('FOURTH_EPOCH',json.dumps(rec),flush=True)
    write('targets.csv',targets);write('references.csv',matches);write('registration.csv',summ)
    if not summ:raise RuntimeError('No F140M registration with >=3 references')
    # Actual JWST barycentric equatorial vectors; no geocentric approximation fallback.
    try:
        from astroquery.jplhorizons import Horizons
        epochs=sorted(set(float(x['mid_mjd']) for x in regs.values())|set(q['mid_mjd'] for q in targets))
        times=Time(epochs,format='mjd',scale='utc')
        vec=Horizons(id='-170',location='500@0',epochs=times.tdb.jd.tolist()).vectors(refplane='earth')
        xyz=np.array([vec[k].quantity.to_value(u.au) for k in ('x','y','z')]).T
        if len(xyz)!=len(epochs) or not np.allclose(np.asarray(vec['datetime_jd']),times.tdb.jd,atol=1e-7,rtol=0):raise RuntimeError('Horizons epoch mismatch')
        ra,dec=np.deg2rad([RA,DEC]);east=np.array([-np.sin(ra),np.cos(ra),0]);north=np.array([-np.cos(ra)*np.sin(dec),-np.sin(ra)*np.sin(dec),np.cos(dec)])
        rows=[dict(mid_mjd=t,px=float(-v@east),py=float(-v@north),x_au=v[0],y_au=v[1],z_au=v[2]) for t,v in zip(epochs,xyz)]
        write('jwst_parallax_factors.csv',rows);print('FACTORS',json.dumps(rows),flush=True)
    except Exception as exc:
        (OUT/'horizons_error.txt').write_text(repr(exc));print('HORIZONS_ERROR',repr(exc),flush=True)
if __name__=='__main__':main()
