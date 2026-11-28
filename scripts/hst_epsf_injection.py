"""HST WFC3/IR moving-source pilot; actual FLT pixels, local registration, visit likelihood.
No validated detection or parallax claim from a selected maximum alone.
"""
import csv,json,math
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote
import requests
import numpy as np
from scipy.special import erf
from scipy.optimize import least_squares
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.mast import Observations
from audit_calibrated_detector_local import exposure_mjd,RA,DEC
OUT=Path('hst_results');OUT.mkdir(exist_ok=True)
DATA=Path('hst_flt');DATA.mkdir(exist_ok=True)
CENTER=SkyCoord(RA,DEC,unit='deg')
VISITS={'ibhm55':('2010-11-03',{'F140W'},['ibhm55030','hst_12328_55_wfc3_ir_f140w_ibhm55']),
'iboa30':('2011-12-10',{'F125W','F160W'},['iboa30030','iboa30040','hst_12440_30_wfc3_ir_f125w_iboa30','hst_12440_30_wfc3_ir_f160w_iboa30']),
'iboa74':('2012-01-26',{'F125W','F160W'},['iboa74030','iboa74bsq','iboa74bwq','hst_12440_74_wfc3_ir_f125w_iboa74','hst_12440_74_wfc3_ir_f160w_iboa74'])}
URL='https://mast.stsci.edu/api/v0.1/Download/file?uri='
def write(name,rows):
    if not rows:return
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=sorted(set().union(*(q.keys() for q in rows))));w.writeheader();w.writerows(rows)
def read(name):
    p=list(Path('baseline').rglob(name))
    if len(p)!=1:raise ValueError((name,p))
    with p[0].open() as f:return list(csv.DictReader(f))
def refs():
    a={q['id']:q for q in read('anchor_refs.csv')}
    regs={q['filename']:q for q in read('registration_pilot.csv')}
    by=defaultdict(lambda:defaultdict(list))
    for q in read('reference_centroids.csv'):
        if q['filename'] not in regs:continue
        r=regs[q['filename']]
        by[q['ref_id']][q['program']].append([float(q['east_mas'])-float(r['east_mas']),float(q['north_mas'])-float(r['north_mas'])])
    result={}
    for k,v in by.items():
        if len(v)!=3:continue
        z=np.array([np.median(x,axis=0) for x in v.values()]);c=np.median(z,axis=0)
        if np.max(np.linalg.norm(z-c,axis=1))>20:continue
        result[k]=SkyCoord(float(a[k]['ra']),float(a[k]['dec']),unit='deg').spherical_offsets_by(c[0]*u.mas,c[1]*u.mas)
    return result
def download(name,uri):
    p=DATA/name
    if p.exists():return p
    with requests.get(URL+quote(uri,safe=':/'),stream=True,timeout=(20,120)) as resp:
        resp.raise_for_status()
        temp=p.with_suffix('.part')
        try:
            with temp.open('wb') as f:
                for c in resp.iter_content(2**20):f.write(c)
            with temp.open('rb') as f:
                if not f.read(20).startswith(b'SIMPLE  ='):raise ValueError('Not FITS')
            temp.rename(p)
        finally:temp.unlink(missing_ok=True)
    return p
def discover(prefix,ids):
    products={};manifest=[]
    for query in [ids,[prefix+'*']]:
        try:
            obs=Observations.query_criteria(obs_id=query)
            for p in Observations.get_product_list(obs):
                name=str(p['productFilename']).lower()
                if name.startswith(prefix) and name.endswith(('_flt.fits','_asn.fits')):
                    products[name]=str(p['dataURI'])
            manifest.append(dict(prefix=prefix,query=str(query),n_obs=len(obs),products=len(products)))
        except Exception as e:manifest.append(dict(prefix=prefix,query=str(query),error=repr(e)))
    for root in ids:
        if root.startswith(prefix) and len(root)==9:
            suffix='_asn.fits' if root.endswith('0') else '_flt.fits'
            name=root+suffix;products.setdefault(name,'mast:HST/product/'+name)
    for name,uri in list(products.items()):
        if not name.endswith('_asn.fits'):continue
        try:
            with fits.open(download(name,uri),memmap=False) as h:
                for q in h[1].data:
                    if 'EXP' not in str(q['MEMTYPE']).upper():continue
                    root=str(q['MEMNAME']).strip().lower()
                    if root.startswith(prefix):
                        n=root+'_flt.fits';products.setdefault(n,'mast:HST/product/'+n)
        except Exception as e:manifest.append(dict(prefix=prefix,query=name,error=repr(e)))
    return {k:v for k,v in products.items() if k.endswith('_flt.fits')},manifest
def offset(a,b):
    x,y=a.spherical_offsets_to(b);return np.array([x.mas,y.mas])
def pixel(w,pos):
    return np.asarray(w.all_world2pix([[pos.ra.deg,pos.dec.deg]],0)[0],float)
def sky(w,x,y):
    return SkyCoord(*w.all_pix2world([[x,y]],0)[0],unit='deg')
def patch(a,e,d,x,y,h):
    ix,iy=round(x),round(y)
    if not(h<=ix<a.shape[1]-h and h<=iy<a.shape[0]-h):raise ValueError('outside safe detector')
    sl=np.s_[iy-h:iy+h+1,ix-h:ix+h+1]
    z=np.asarray(a[sl],float);err=np.asarray(e[sl],float);dq=d[sl]
    yy,xx=np.indices(z.shape,dtype=float);xx+=ix-h;yy+=iy-h
    good=np.isfinite(z)&np.isfinite(err)&(err>0)&(dq==0)
    if good.mean()<.65:raise ValueError('too many masked pixels')
    skyann=good&(np.hypot(xx-x,yy-y)>h*.65)
    if skyann.sum()<20:raise ValueError('insufficient sky')
    bg=np.median(z[skyann]);noise=1.4826*np.median(abs(z[skyann]-bg))
    return z, np.maximum(err,max(noise,1e-10)),good,xx,yy,bg
def template(xx,yy,x,y,s):
    def axis(z,t):return .5*(erf((z-t+.5)/(np.sqrt(2)*s))-erf((z-t-.5)/(np.sqrt(2)*s)))
    return axis(xx,x)*axis(yy,y)
def ref_fit(a,e,d,w,pos):
    x,y=pixel(w,pos);z,err,g,xx,yy,bg=patch(a,e,d,x,y,7)
    def res(p):
        amp,dx,dy,s,c,gx,gy=p
        return ((amp*template(xx,yy,x+dx,y+dy,s)+c+gx*(xx-x)+gy*(yy-y)-z)/err)[g]
    amp=max(np.sum((z-bg)[g&(np.hypot(xx-x,yy-y)<2)]),np.median(err[g])*10)
    lo=[0,-3,-3,.35,-np.inf,-np.inf,-np.inf];hi=[np.inf,3,3,2.,np.inf,np.inf,np.inf]
    rr=least_squares(res,[amp,0,0,.65,bg,0,0],bounds=(lo,hi),max_nfev=150)
    p=rr.x;cov=np.linalg.pinv(rr.jac.T@rr.jac)*max(1,2*rr.cost/(g.sum()-7))
    snr=p[0]/max(np.sqrt(max(cov[0,0],0)),1e-20)
    if not rr.success or np.any(rr.active_mask) or snr<8 or np.sqrt(np.max(np.diag(cov)[1:3]))>.35:raise ValueError('low quality reference fit')
    return offset(pos,sky(w,x+p[1],y+p[2])),float(snr)
def pred(t):
    # Based on the four measured NIRCam positions; approximate search prior only.
    dt=(t-60308.94)/365.25
    return CENTER.spherical_offsets_by((2.57-70.53*dt)*u.mas,(17.59-61.06*dt)*u.mas)
def scan(a,e,d,w,position,grid,sigma):
    x,y=pixel(w,position);z,err,g,xx,yy,bg=patch(a,e,d,x,y,12)
    # Project constant and linear background out of whitened data/template.
    D=np.stack([np.ones(z.shape),xx-x,yy-y],-1)[g]/err[g,None]
    Q=np.linalg.qr(D,mode='reduced')[0];v=z[g]/err[g];v-=Q@(Q.T@v)
    pts=position.spherical_offsets_by(grid[:,0]*u.arcsec,grid[:,1]*u.arcsec)
    xy=w.all_world2pix(np.column_stack([pts.ra.deg,pts.dec.deg]),0)
    scores=np.full(len(grid),np.nan);flux=scores.copy();variance=scores.copy()
    for i,(px,py) in enumerate(xy):
        if not np.isfinite(px+py):continue
        if np.sum(g&(np.hypot(xx-px,yy-py)<1.8))<7:continue
        k=template(xx,yy,px,py,sigma)[g]/err[g];k-=Q@(Q.T@k)
        norm=k@k
        if norm<=0:continue
        scores[i]=k@v;flux[i]=scores[i]/norm;variance[i]=1/norm
    return flux,variance
def combined(group):
    f=np.array([q[0] for q in group]);v=np.array([q[1] for q in group])
    ok=np.isfinite(f)&np.isfinite(v)&(v>0)
    weights=np.where(ok,1/v,0);numerator=np.sum(np.where(ok,f*weights,0),axis=0);den=weights.sum(axis=0)
    sn=np.divide(numerator,np.sqrt(den),out=np.full_like(den,np.nan),where=den>0)
    sn[ok.sum(axis=0)!=len(group)]=np.nan
    return sn

def injection_audit(a,e,d,w,center,sigma,amp,filt,name):
    from photutils.psf import GriddedPSFModel
    from astropy.nddata import NDData
    from scipy.optimize import least_squares
    psf=GriddedPSFModel.read('psf_'+filt+'.fits',format='stdpsf')
    x,y=pixel(w,center)
    z,err,g,xx,yy,bg=patch(a,e,d,x,y,8)
    def resid(p):
        f,dx,dy,b,gx,gy=p
        psf.x_0=x+dx;psf.y_0=y+dy;psf.flux=f
        return ((psf(xx,yy)+b+gx*(xx-x)+gy*(yy-y)-z)/err)[g]
    fit=least_squares(resid,[amp,0,0,bg,0,0],bounds=([0,-5,-5,-np.inf,-np.inf,-np.inf],[np.inf,5,5,np.inf,np.inf,np.inf]))
    v=offset(center,sky(w,x+fit.x[1],y+fit.x[2]))
    print('TARGET_EPSF',json.dumps(dict(filename=name,filter=filt,east=float(v[0]),north=float(v[1]),flux=float(fit.x[0]),valid=bool(fit.success and not np.any(fit.active_mask)),chi2=float(2*fit.cost),dof=int(g.sum()-6))),flush=True)
    rng=np.random.default_rng(639188)
    rows=[]
    for j in range(48):
        ang=2*np.pi*j/48
        pos=center.spherical_offsets_by((3*np.cos(ang))*u.arcsec,(3*np.sin(ang))*u.arcsec)
        x,y=pixel(w,pos);x+=rng.uniform(-.5,.5);y+=rng.uniform(-.5,.5)
        try:
            z,err,g,xx,yy,bg=patch(a,e,d,x,y,6)
            psf.x_0=x;psf.y_0=y;psf.flux=amp
            inj=z+psf(xx,yy)
            def fit(model):
                def res(p):
                    f,dx,dy,b,gx,gy=p
                    if model=='epsf':
                        psf.x_0=x+dx;psf.y_0=y+dy;psf.flux=f;k=psf(xx,yy)
                    else:k=f*template(xx,yy,x+dx,y+dy,sigma)
                    return ((k+b+gx*(xx-x)+gy*(yy-y)-inj)/err)[g]
                r=least_squares(res,[amp,0,0,bg,0,0],bounds=([0,-1,-1,-np.inf,-np.inf,-np.inf],[np.inf,1,1,np.inf,np.inf,np.inf]))
                v=offset(sky(w,x,y),sky(w,x+r.x[1],y+r.x[2]))
                return dict(east=float(v[0]),north=float(v[1]),valid=bool(r.success and not np.any(r.active_mask)),flux=float(r.x[0]))
            rows.append(dict(filename=name,trial=j,filter=filt,gaussian=fit('gaussian'),epsf=fit('epsf')))
        except Exception as ex:rows.append(dict(filename=name,trial=j,error=repr(ex)))
    (OUT/(name+'_injection.json')).write_text(json.dumps(rows))
    for model in ['gaussian','epsf']:
        z=np.array([[r[model]['east'],r[model]['north']] for r in rows if model in r and r[model]['valid']])
        if len(z)>2:print('INJECTION',json.dumps(dict(filename=name,model=model,n=len(z),bias=z.mean(0).tolist(),cov=np.cov(z.T).tolist())),flush=True)

def main():
    ref=refs();print('NIRCAM_REFERENCE_CANDIDATES',len(ref),flush=True)
    audit=[];refrows=[];queries=[];groups=defaultdict(list);blank=defaultdict(lambda:defaultdict(list))
    xy=np.arange(-.60,.601,.04);xx,yy=np.meshgrid(xy,xy);grid=np.column_stack([xx.ravel(),yy.ravel()])
    blank_offsets=[(3*np.cos(p),3*np.sin(p)) for p in np.linspace(0,2*np.pi,12,endpoint=False)]
    for prefix,(date,filters,ids) in VISITS.items():
        products,m=discover(prefix,ids);queries+=m;write('query_manifest.csv',queries)
        for name,uri in sorted(products.items()):
            row=dict(filename=name,visit=date,status='error')
            try:
                path=download(name,uri)
                with fits.open(path,memmap=False) as h:
                    ph=h[0].header;filt=str(ph.get('FILTER','')).upper();row['filter']=filt
                    if filt not in filters or ph.get('DETECTOR')!='IR':row['status']='other_filter_or_detector';continue
                    sh=h['SCI',1];a=sh.data;e=h['ERR',1].data;d=h['DQ',1].data
                    t=exposure_mjd(ph,sh.header);w=WCS(sh.header,h)
                    position=pred(t);x,y=pixel(w,position)
                    row.update(mid_mjd=t,pred_ra=position.ra.deg,pred_dec=position.dec.deg,pred_x=x,pred_y=y,bunit=sh.header.get('BUNIT',''),exptime=ph.get('EXPTIME'))
                    if not (12<x<a.shape[1]-13 and 12<y<a.shape[0]-13):row['status']='outside_detector';continue
                    offsets=[]
                    for k,p in ref.items():
                        try:
                            v,s=ref_fit(a,e,d,w,p);offsets.append(v)
                            refrows.append(dict(filename=name,visit=date,filter=filt,ref_id=k,east_mas=v[0],north_mas=v[1],snr=s))
                        except Exception:continue
                    z=np.array(offsets)
                    if len(z)>=4:
                        c=np.median(z,axis=0);mad=np.maximum(15,1.4826*np.median(abs(z-c),axis=0));z=z[np.all(abs(z-c)<=4*mad,axis=1)]
                    if len(z)>=4:
                        c=np.median(z,axis=0);loo=np.array([np.median(np.delete(z,i,axis=0),axis=0) for i in range(len(z))])
                        row.update(n_refs=len(z),translation_east_mas=c[0],translation_north_mas=c[1],ref_mad_mas=(1.4826*np.median(abs(z-c),axis=0)).tolist(),loo_range_mas=np.ptp(loo,axis=0).tolist(),registered=True)
                    else:
                        c=np.zeros(2);row.update(n_refs=len(z),registered=False)
                    center=position.spherical_offsets_by(c[0]*u.mas,c[1]*u.mas)
                    sigma={'F125W':.55,'F140W':.60,'F160W':.65}[filt]
                    ff,vv=scan(a,e,d,w,center,grid,sigma);sn=combined([(ff,vv)])
                    i=np.nanargmax(sn);row.update(status='pixel_likelihood',max_snr=float(sn[i]),dx_arcsec=grid[i,0],dy_arcsec=grid[i,1],flux_native=float(ff[i]),flux_sigma_native=float(np.sqrt(vv[i])),boundary=bool(np.max(abs(grid[i]))>.56))
                    injection_audit(a,e,d,w,center.spherical_offsets_by(grid[i,0]*u.arcsec,grid[i,1]*u.arcsec),sigma,float(ff[i]),filt,name)
                    key=(date,filt,bool(row['registered']))
                    np.savez_compressed(OUT/f'{name}_likelihood.npz',grid=grid,flux=ff,variance=vv)
                    groups[key].append((ff,vv))
                    for j,(bx,by) in enumerate(blank_offsets):
                        try:
                            bpos=center.spherical_offsets_by(bx*u.arcsec,by*u.arcsec)
                            blank[key][j].append(scan(a,e,d,w,bpos,grid,sigma))
                        except Exception:continue
            except Exception as ex:row['error']=repr(ex)
            finally:
                audit.append(row);write('exposure_audit.csv',audit);write('references.csv',refrows);print('HST_FLT',json.dumps(row,default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist()),flush=True)
    summary=[]
    for key,g in groups.items():
        sn=combined(g);i=np.nanargmax(sn);controls=[]
        for j,b in blank[key].items():
            if len(b)==len(g):
                s=combined(b)
                if np.isfinite(s).any():controls.append(float(np.nanmax(s)))
        date,filt,reg=key
        r=dict(visit=date,filter=filt,registered=reg,n_exposures=len(g),peak_snr=float(sn[i]),dx_arcsec=grid[i,0],dy_arcsec=grid[i,1],boundary=bool(np.max(abs(grid[i]))>.56),n_control_regions=len(controls),control_max_snr=controls,controls_ge_peak=int(sum(bool(x>=sn[i]) for x in controls)),stationary_verified=False,source_association_validated=False)
        summary.append(r);print('HST_VISIT',json.dumps(r,default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist()),flush=True)
        np.savez_compressed(OUT/f'{date}_{filt}_{reg}_likelihood.npz',grid=grid,snr=sn,controls=controls)
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist()))
    if not summary:raise RuntimeError('No HST individual detector likelihood computed; inspect audit')
if __name__=='__main__':main()
