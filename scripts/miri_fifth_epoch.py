"""Fifth epoch pilot: reuse actual MIRI target fits, tie to NIRCam common controls.
Cross-filter/instrument offsets are diagnostics; never mark an epoch validated here.
"""
import csv,json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
from scipy.special import erf
from scipy.optimize import least_squares
from jwst import datamodels
from nircam_same_filter_astrometry import download,offset
from audit_calibrated_detector_local import exposure_mjd,RA,DEC
OUT=Path('fifth_epoch');OUT.mkdir(exist_ok=True)
CENTER=SkyCoord(RA,DEC,unit='deg')

def read(root,name):
    p=list(Path(root).rglob(name))
    if len(p)!=1:raise RuntimeError(f'Expected exactly one {name}: {p}')
    with p[0].open() as f:return list(csv.DictReader(f))
def write(name,rows):
    if rows:
        with (OUT/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=sorted(set().union(*(q.keys() for q in rows))));w.writeheader();w.writerows(rows)
def baseline():
    anchors={q['id']:q for q in read('baseline','anchor_refs.csv')}
    regs={q['filename']:q for q in read('baseline','registration_pilot.csv')}
    by={}
    for q in read('baseline','reference_centroids.csv'):
        if q['filename'] not in regs:continue
        r=regs[q['filename']]
        by.setdefault(q['ref_id'],{}).setdefault(q['program'],[]).append([float(q['east_mas'])-float(r['east_mas']),float(q['north_mas'])-float(r['north_mas'])])
    result={}
    for k,v in by.items():
        if len(v)!=3:continue
        a=np.array([np.median(x,axis=0) for x in v.values()])
        c=np.median(a,axis=0)
        pos=SkyCoord(float(anchors[k]['ra']),float(anchors[k]['dec']),unit='deg').spherical_offsets_by(c[0]*u.mas,c[1]*u.mas)
        result[k]=(pos,float(np.max(np.linalg.norm(a-c,axis=1))))
    return result

def control(a,e,d,w,pos,filt):
    x,y=map(float,w.invert(pos.ra.deg,pos.dec.deg))
    if not np.isfinite(x+y) or not (12<x<a.shape[1]-13 and 12<y<a.shape[0]-13):raise ValueError('off detector')
    ix,iy=round(x),round(y);h=10
    sl=np.s_[iy-h:iy+h+1,ix-h:ix+h+1]
    z=np.asarray(a[sl],float);err=np.asarray(e[sl],float);dq=d[sl]
    yy,xx=np.indices(z.shape,dtype=float);xx+=ix-h-x;yy+=iy-h-y
    rad=np.hypot(xx,yy);good=np.isfinite(z)&np.isfinite(err)&(err>0)&((dq&3)==0)
    core=rad<4
    if good[core].mean()<.9 or good.sum()<300:raise ValueError('masked')
    sky=good&(rad>7);bg=np.median(z[sky]);floor=max(1e-12,1.4826*np.median(abs(z[sky]-bg)))
    err=np.maximum(err,floor)
    def model(p):
        amp,dx,dy,sx,sy,c,g,h=p
        def axis(z,t,s):return .5*(erf((z-t+.5)/(np.sqrt(2)*s))-erf((z-t-.5)/(np.sqrt(2)*s)))
        return amp*axis(xx,dx,sx)*axis(yy,dy,sy)+c+g*xx+h*yy
    def res(p):return ((model(p)-z)/err)[good]
    sigma={'F770W':1.05,'F1000W':1.25,'F1280W':1.65,'F1500W':1.9}[filt]
    lo=[0,-2.5,-2.5,.65,.65,-np.inf,-np.inf,-np.inf];hi=[np.inf,2.5,2.5,4.,4.,np.inf,np.inf,np.inf]
    amp=max(floor*10,np.sum(z[good&core]-bg))
    fit=least_squares(res,[amp,0,0,sigma,sigma,bg,0,0],bounds=(lo,hi),max_nfev=200)
    p=fit.x
    X=np.stack([np.ones_like(xx),xx,yy],axis=-1)[good]/err[good,None]
    v=z[good]/err[good];null=v-X@np.linalg.lstsq(X,v,rcond=None)[0]
    delta=float(np.sum(null**2)-2*fit.cost)
    cov=np.linalg.pinv(fit.jac.T@fit.jac)*max(1,2*fit.cost/(good.sum()-8))
    formal=np.sqrt(np.maximum(0,np.diag(cov)[1:3]))
    if not fit.success or np.any(fit.active_mask) or delta<50 or np.max(formal)>.5:raise ValueError('control fit quality')
    ps=SkyCoord(*map(float,w(x+p[1],y+p[2])),unit='deg')
    dx,dy=offset(pos,ps)
    return dict(east=dx,north=dy,delta_chi2=delta,formal_x_pix=formal[0],formal_y_pix=formal[1],sigma_x_pix=p[3],sigma_y_pix=p[4],psf_approx='integrated_gaussian')

def main():
    base=baseline();print('BASELINE_COMMON',len(base),flush=True)
    diag={q['filename']:q for q in read('diagnostic','miri_cal_detector_diagnostic.csv')}
    psf=read('psf','miri_stpsf_diagnostic.csv')
    allrefs=[];targets=[];regs=[];errors=[]
    for q in psf:
        f=q['filter']
        if f not in ('F1280W','F1500W'):continue
        name=q['filename']
        try:
            if q['free_optimizer_success'].lower()!='true' or q['free_near_search_boundary'].lower()=='true':raise ValueError('target fit quality')
            d=diag[name]
            x=float(d['x_photometry'])+float(q['free_dx_pix']);y=float(d['y_photometry'])+float(q['free_dy_pix'])
            fn=download('mast:JWST/product/'+name,name)
            with fits.open(fn,memmap=False) as hd,datamodels.open(fn) as m:
                if hd[0].header['FILTER']!=f:raise ValueError('filter mismatch')
                a,e,dq=hd['SCI'].data,hd['ERR'].data,hd['DQ'].data;w=m.meta.wcs
                t=exposure_mjd(hd[0].header,hd['SCI'].header)
                p=SkyCoord(*map(float,w(x,y)),unit='deg');raw=np.array(offset(CENTER,p))
                targets.append(dict(filename=name,filter=f,mid_mjd=t,raw_east_mas=raw[0],raw_north_mas=raw[1],flux_ujy=float(q['free_flux_ujy']),x=x,y=y))
                for k,(pos,scatter) in base.items():
                    try:
                        z=control(a,e,dq,w,pos,f)
                        allrefs.append(dict(filename=name,filter=f,ref_id=k,baseline_scatter_mas=scatter,**z))
                    except Exception as ex:errors.append(dict(filename=name,ref_id=k,error=str(ex)))
                for cut in (10,20,50):
                    rr=[z for z in allrefs if z['filename']==name and z['baseline_scatter_mas']<=cut]
                    if len(rr)<3:continue
                    v=np.array([[z['east'],z['north']] for z in rr]);c=np.median(v,axis=0)
                    mad=np.maximum(1,1.4826*np.median(abs(v-c),axis=0));keep=np.all(abs(v-c)<=4*mad,axis=1)
                    if keep.sum()<3:continue
                    v=v[keep];c=np.median(v,axis=0);loo=np.array([np.median(np.delete(v,i,axis=0),axis=0) for i in range(len(v))])
                    rec=dict(filename=name,filter=f,mid_mjd=t,n_ref=len(v),cut_mas=cut,registered_east_mas=float(raw[0]-c[0]),registered_north_mas=float(raw[1]-c[1]),translation_east_mas=float(c[0]),translation_north_mas=float(c[1]),reference_mad_mas=(1.4826*np.median(abs(v-c),axis=0)).tolist(),loo_range_mas=np.ptp(loo,axis=0).tolist(),stationary_verified=False,cross_filter_validated=False)
                    regs.append(rec);print('FIFTH_REG',json.dumps(rec),flush=True)
        except Exception as ex:errors.append(dict(filename=name,error=repr(ex)))
        write('targets.csv',targets);write('references.csv',allrefs);write('registration.csv',regs);write('errors.csv',errors)
    print('FIFTH_COUNTS',json.dumps(dict(targets=len(targets),reference_rows=len(allrefs),registered_files=len(set(q['filename'] for q in regs)))),flush=True)
    if not targets:raise RuntimeError('No usable 2025 target data')
    from astroquery.jplhorizons import Horizons
    dates=sorted(set(q['mid_mjd'] for q in targets))
    times=Time(dates,format='mjd',scale='utc')
    vec=Horizons(id='-170',location='500@0',epochs=times.tdb.jd.tolist()).vectors(refplane='earth')
    xyz=np.array([vec[k].quantity.to_value(u.au) for k in ('x','y','z')]).T
    if len(xyz)!=len(dates) or not np.allclose(np.asarray(vec['datetime_jd']),times.tdb.jd,atol=1e-7,rtol=0):raise RuntimeError('Horizons date mismatch')
    ra,dec=np.deg2rad([RA,DEC]);east=np.array([-np.sin(ra),np.cos(ra),0]);north=np.array([-np.cos(ra)*np.sin(dec),-np.sin(ra)*np.sin(dec),np.cos(dec)])
    fac=[dict(mid_mjd=t,px=float(-v@east),py=float(-v@north),x_au=v[0],y_au=v[1],z_au=v[2]) for t,v in zip(dates,xyz)]
    write('factors.csv',fac);print('FIFTH_FACTORS',json.dumps(fac),flush=True)
    # Equal exposure counts are not independent epochs. Compare filter medians and per-source offsets.
    summary=[]
    for f in ('F1280W','F1500W'):
        for cut in (10,20,50):
            v=[q for q in regs if q['filter']==f and q['cut_mas']==cut]
            if v:
                a=np.array([[q['registered_east_mas'],q['registered_north_mas']] for q in v])
                summary.append(dict(filter=f,cut_mas=cut,n=len(v),median_mas=np.median(a,axis=0).tolist(),sample_sd_mas=np.std(a,axis=0,ddof=1).tolist() if len(v)>1 else None))
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2));print('FIFTH_SUMMARY',json.dumps(summary),flush=True)
    if not regs:raise RuntimeError('No cross-instrument registration; inspect errors and target rows')
if __name__=='__main__':main()
