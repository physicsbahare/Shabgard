"""Same-filter NIRCam detector astrometry pilot; no parallax fit or verified epoch claim."""
import csv,json,math,signal
from pathlib import Path
from urllib.parse import quote
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.mast import Observations
from scipy.ndimage import maximum_filter, gaussian_filter
from jwst import datamodels
import requests
from localize_f140m import fit
from audit_calibrated_detector_local import exposure_mjd,sky_at,RA,DEC

OUT=Path('nircam_results'); OUT.mkdir(exist_ok=True)
DATA=Path('nircam_cal'); DATA.mkdir(exist_ok=True)
ASSOC={'1837':'jw01837-o003_t003_nircam_clear-f115w','1727':'jw01727-o051_t051_nircam_clear-f115w','5893':'jw05893-o014_t014_nircam_clear-f115w'}
URL='https://mast.stsci.edu/api/v0.1/Download/file?uri='
def dump(name,rows):
    if not rows:return
    with (OUT/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=sorted(set().union(*(x.keys() for x in rows))))
        w.writeheader();w.writerows(rows)
def timeout(s,f): raise TimeoutError('remote header timeout')
signal.signal(signal.SIGALRM,timeout)
def get_header(uri):
    signal.alarm(35)
    try:
        with fits.open(URL+quote(uri,safe=':/'),use_fsspec=True,fsspec_kwargs={'block_size':65536},lazy_load_hdus=True) as h:
            return h[0].header.copy(),h['SCI'].header.copy()
    finally:signal.alarm(0)
def footprint(p,s):
    # Header S_REGION polygon only for preselection; actual gWCS is definitive.
    from matplotlib.path import Path as Polygon
    a=str(s.get('S_REGION',p.get('S_REGION',''))).split()
    if len(a)<8 or a[0].upper()!='POLYGON':return None
    try: v=np.array(list(map(float,a[2:]))).reshape(-1,2)
    except ValueError:return None
    v[:,0]=(v[:,0]-RA)*3600*math.cos(math.radians(DEC));v[:,1]=(v[:,1]-DEC)*3600
    return bool(Polygon(v).contains_point((0,0),radius=6.))
def download(uri,filename):
    p=DATA/filename
    if p.exists():return p
    tmp=p.with_suffix('.part')
    try:
        with requests.get(URL+quote(uri,safe=':/'),stream=True,timeout=(30,120)) as resp:
            resp.raise_for_status()
            with tmp.open('wb') as f:
                for chunk in resp.iter_content(2**20):f.write(chunk)
        tmp.rename(p)
    finally:
        if tmp.exists():tmp.unlink()
    return p
def usable(q,threshold=25):
    return q['converged'] and not q['boundary'] and q['delta_chi2']>=threshold and q['formal_x']<.5 and q['formal_y']<.5
def coord(w,q):
    ra,dec=w(q['x'],q['y']);return SkyCoord(float(ra),float(dec),unit='deg')
def offset(ref,pos):
    dx,dy=ref.spherical_offsets_to(pos);return float(dx.mas),float(dy.mas)
def fit_one(a,e,d,w,ra,dec):
    x,y=map(float,w.invert(ra,dec))
    if not(np.isfinite(x+y) and 9<x<a.shape[1]-10 and 9<y<a.shape[0]-10):raise ValueError('outside detector')
    q=fit(a,e,d,x,y);return q,coord(w,q),x,y
def refs_from_anchor(a,e,d,w,tx,ty):
    # Sources in a bounded local field, not a catalog of verified stationary objects.
    h=650; ix,iy=round(tx),round(ty);x0=max(9,ix-h);x1=min(a.shape[1]-9,ix+h);y0=max(9,iy-h);y1=min(a.shape[0]-9,iy+h)
    z=np.asarray(a[y0:y1,x0:x1],float);valid=np.isfinite(z)&((d[y0:y1,x0:x1]&3)==0)
    bg=np.median(z[valid]);noise=max(1e-12,1.4826*np.median(abs(z[valid]-bg)))
    zz=gaussian_filter(np.where(valid,z-bg,0.),.8)
    yy,xx=np.nonzero((zz==maximum_filter(zz,size=9))&(zz>8*noise)&valid)
    order=np.argsort(zz[yy,xx])[::-1];refs=[]
    for i in order[:300]:
        x,y=float(xx[i]+x0),float(yy[i]+y0)
        if np.hypot(x-tx,y-ty)<25 or any(np.hypot(x-q['x'],y-q['y'])<25 for q in refs):continue
        try:
            q=fit(a,e,d,x,y)
            if not usable(q,100):continue
            pos=coord(w,q);refs.append(dict(id=len(refs),ra=pos.ra.deg,dec=pos.dec.deg,x=q['x'],y=q['y']))
        except Exception:continue
        if len(refs)>=70:break
    return refs
def main():
    manifest=[]; files=[]; targets=[]; reference_rows=[];refs=[]
    # Query exact association IDs; recover members from association JSON, with product list fallback.
    for program,assoc in ASSOC.items():
        try:
            obs=Observations.query_criteria(obs_id=assoc,obs_collection='JWST')
            if len(obs)==0:raise RuntimeError('exact association lookup empty')
            products=Observations.get_product_list(obs)
            members=set()
            for row in products:
                name=str(row['productFilename']);uri=str(row['dataURI'])
                if not name.endswith('_asn.json'):continue
                try:
                    resp=requests.get(URL+quote(uri,safe=':/'),timeout=45);resp.raise_for_status()
                    for prod in resp.json().get('products',[]):
                        if prod.get('name')!=assoc:continue
                        for m in prod.get('members',[]):
                            name=Path(m['expname']).name
                            if m.get('exptype')=='science' and '_nrc' in name:
                                members.add(name.replace('_rate.fits','_cal.fits'))
                except Exception as exc:print('ASN error',name,repr(exc),flush=True)
            if not members:
                members={str(q['productFilename']) for q in products if str(q['productFilename']).endswith('_cal.fits') and '_nrc' in str(q['productFilename'])}
            if not members:raise RuntimeError('no CAL members found')
            print('MEMBERS',program,len(members),flush=True)
            # Limit downloads only after header footprint screening; flag any unprocessed members.
            accepted=0
            for name in sorted(members):
                if not name.endswith('_cal.fits'):continue
                uri='mast:JWST/product/'+name
                rec=dict(program=program,association=assoc,filename=name,status='error')
                try:
                    p,s=get_header(uri)
                    if p.get('FILTER')!='F115W':rec['status']='other_filter';continue
                    hit=footprint(p,s);rec['header_footprint']=str(hit)
                    if hit is False:rec['status']='outside_header_polygon';continue
                    if accepted>=24:rec['status']='download_cap_not_tested';continue
                    fn=download(uri,name);accepted+=1;files.append((program,fn));rec['status']='downloaded'
                except Exception as exc:rec['error']=repr(exc)
                finally:manifest.append(rec);dump('member_audit.csv',manifest)
        except Exception as exc:
            manifest.append(dict(program=program,association=assoc,status='association_error',error=repr(exc)));dump('member_audit.csv',manifest)
    for program,fn in files:
        row=dict(program=program,filename=fn.name,status='error')
        try:
            with fits.open(fn,memmap=False) as h,datamodels.open(fn) as m:
                p,s=h[0].header,h['SCI'].header;a,e,d=h['SCI'].data,h['ERR'].data,h['DQ'].data;w=m.meta.wcs
                t=exposure_mjd(p,s);row['mid_mjd']=t
                for key in ('MJD-AVG','EXPSTART','EXPEND','JWST_X','JWST_Y','JWST_Z','OBS_VX','OBS_VY','OBS_VZ','EPH_TIME','REFFRAME'):
                    if key in p:row[key]=p[key]
                # Provisional catalog epoch used solely for locating target, not parallax subtraction.
                ra,dec=sky_at(t,60308)
                q,pos,x,y=fit_one(a,e,d,w,ra,dec)
                row.update(q);row.update(ra_deg=pos.ra.deg,dec_deg=pos.dec.deg,status='candidate' if usable(q) else 'low_quality',pred_x=x,pred_y=y)
                if usable(q) and not refs:
                    refs=refs_from_anchor(a,e,d,w,x,y);dump('anchor_refs.csv',refs)
                    print('ANCHOR',fn.name,len(refs),flush=True)
                for ref in refs:
                    try:
                        z,ps,_,_=fit_one(a,e,d,w,ref['ra'],ref['dec'])
                        if not usable(z,50):continue
                        dx,dy=offset(SkyCoord(ref['ra'],ref['dec'],unit='deg'),ps)
                        reference_rows.append(dict(program=program,filename=fn.name,ref_id=ref['id'],east_mas=dx,north_mas=dy,**z))
                    except Exception:continue
        except Exception as exc:row['error']=repr(exc)
        targets.append(row);dump('target_centroids.csv',targets);dump('reference_centroids.csv',reference_rows)
        print('TARGET',json.dumps(row),flush=True)
    summary=[]
    for row in targets:
        if row['status']!='candidate':continue
        a=np.array([[z['east_mas'],z['north_mas']] for z in reference_rows if z['filename']==row['filename']])
        if len(a)<6:continue
        center=np.median(a,axis=0)
        mad=np.maximum(1.,1.4826*np.median(abs(a-center),axis=0))
        good=np.all(abs(a-center)<=4*mad,axis=1);b=a[good]
        if len(b)<6:continue
        c=np.median(b,axis=0);scatter=1.4826*np.median(abs(b-c),axis=0)
        pos=SkyCoord(row['ra_deg'],row['dec_deg'],unit='deg').spherical_offsets_by(-c[0]*u.mas,-c[1]*u.mas)
        summary.append(dict(filename=row['filename'],program=row['program'],mid_mjd=row['mid_mjd'],n_ref=len(b),east_mas=c[0],north_mas=c[1],scatter_east_mas=scatter[0],scatter_north_mas=scatter[1],registered_ra_deg=pos.ra.deg,registered_dec_deg=pos.dec.deg,stationary_verified=False))
    dump('registration_pilot.csv',summary)
    print('REGISTRATION',json.dumps(summary),flush=True)
    if not summary:raise RuntimeError('no diagnostic registrations')
if __name__=='__main__':main()
