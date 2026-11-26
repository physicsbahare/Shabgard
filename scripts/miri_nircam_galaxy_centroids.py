"""Measure candidate galaxy aperture centroids directly in F140M and F770W CAL.
Signed flux moments, local plane, matched angular apertures, no stationary classification.
"""
import argparse,csv,glob,json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from jwst import datamodels
from audit_calibrated_detector_local import exposure_mjd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

IDS={'638656','638807','639322','639535','639727','639977'}
def measure(a,e,dq,x,y,jac,r):
    # Local tangent plane Jacobian in arcsec/pixel; use common sky radii.
    half=int(np.ceil((r+0.55)/np.linalg.svd(jac,compute_uv=False).min()))+2
    ix,iy=int(round(x)),int(round(y))
    if ix-half<0 or iy-half<0 or ix+half>=a.shape[1] or iy+half>=a.shape[0]: raise ValueError('edge')
    sl=np.s_[iy-half:iy+half+1,ix-half:ix+half+1]
    z=np.asarray(a[sl],float);err=np.asarray(e[sl],float);q=np.asarray(dq[sl])
    yy,xx=np.indices(z.shape,dtype=float);xx+=ix-half-x;yy+=iy-half-y
    dx=jac[0,0]*xx+jac[0,1]*yy;dy=jac[1,0]*xx+jac[1,1]*yy
    rad=np.hypot(dx,dy);good=np.isfinite(z)&np.isfinite(err)&(err>0)&((q&3)==0)
    ap=rad<=r;ann=(rad>=r+0.15)&(rad<=r+0.5)
    if good[ap].mean()<.95 or good[ann].mean()<.8: raise ValueError('masked aperture or sky')
    aa=ap&good;bb=ann&good
    X=np.stack([np.ones_like(z),dx,dy],axis=-1)
    M=X[bb];v=z[bb];b=np.linalg.lstsq(M,v,rcond=None)[0]
    for _ in range(4):
        res=v-M@b;med=np.median(res);s=max(1e-12,1.4826*np.median(abs(res-med)))
        keep=abs(res-med)<=3*s
        if keep.sum()<20: break
        b=np.linalg.lstsq(M[keep],v[keep],rcond=None)[0]
    v=z-X@b;fl=v[aa].sum();noise=np.sqrt(np.sum(err[aa]**2))
    if not np.isfinite(fl) or fl<=0: raise ValueError('nonpositive aperture')
    cen=np.array([(v[aa]*dx[aa]).sum(),(v[aa]*dy[aa]).sum()])/fl
    if np.linalg.norm(cen)>r: raise ValueError('centroid outside aperture')
    return cen,float(fl/noise),v,good,half

def main():
    p=argparse.ArgumentParser();p.add_argument('--refs',required=True);p.add_argument('--cal',required=True);p.add_argument('--out',required=True);args=p.parse_args()
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    refs=[r for r in csv.DictReader(open(args.refs)) if r['catalog_id'] in IDS and r['rank']=='1' and float(r['separation_arcsec'])<.3]
    if len(refs)!=6 or len({r['catalog_id'] for r in refs})!=6: raise ValueError('Expected six vetted candidate IDs')
    rows=[]
    for fn in sorted(glob.glob(str(Path(args.cal)/'*_cal.fits'))):
        with fits.open(fn,memmap=False) as h,datamodels.open(fn) as m:
            a,e,q=h['SCI'].data,h['ERR'].data,h['DQ'].data;w=m.meta.wcs
            filt=h[0].header['FILTER'];t=exposure_mjd(h[0].header,h['SCI'].header)
            fig,axs=plt.subplots(2,3,figsize=(10,7));nplot=0
            for ax,ref in zip(axs.flat,refs):
                ax.set_title(ref['catalog_id']);ax.set_axis_off()
                pos=SkyCoord(float(ref['catalog_ra_deg']),float(ref['catalog_dec_deg']),unit='deg')
                x,y=map(float,w.invert(pos.ra.deg,pos.dec.deg))
                if not np.isfinite(x+y): continue
                jac=np.zeros((2,2))
                for j in range(2):
                    d=np.zeros(2);d[j]=.05
                    s1=SkyCoord(*map(float,w(x+d[0],y+d[1])),unit='deg');s0=SkyCoord(*map(float,w(x-d[0],y-d[1])),unit='deg')
                    dx,dy=s0.spherical_offsets_to(s1);jac[:,j]=[dx.arcsec/.1,dy.arcsec/.1]
                for r in [.35,.5,.7]:
                    rec=dict(filename=Path(fn).name,filter=filt,mid_mjd=t,catalog_id=ref['catalog_id'],source_id=ref['source_id'],radius_arcsec=r,status='error',stationary_verified=False)
                    try:
                        cen,snr,v,good,half=measure(a,e,q,x,y,jac,r)
                        sky=pos.spherical_offsets_by(cen[0]*u.arcsec,cen[1]*u.arcsec)
                        rec.update(status='diagnostic',east_mas=cen[0]*1000,north_mas=cen[1]*1000,ra_deg=sky.ra.deg,dec_deg=sky.dec.deg,snr_proxy=snr)
                        if r==.5:
                            z=np.where(good,v,np.nan);lo,hi=np.nanpercentile(z,[5,99]);ax.imshow(z,origin='lower',vmin=lo,vmax=hi,cmap='gray');ax.plot(half+x-round(x),half+y-round(y),'r+');nplot+=1
                    except Exception as ex: rec['error']=str(ex)
                    rows.append(rec)
            if nplot: fig.tight_layout();fig.savefig(out/(Path(fn).stem+'_refs.png'),dpi=110)
            plt.close(fig)
    if not rows: raise ValueError('No candidate measurements')
    with (out/'galaxy_aperture_centroids.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=sorted(set().union(*(r.keys() for r in rows))));w.writeheader();w.writerows(rows)
    summary={}
    for filt in sorted({r['filter'] for r in rows}):
        r=[q for q in rows if q['filter']==filt and q['status']=='diagnostic' and q['radius_arcsec']==.5 and q['snr_proxy']>5]
        summary[filt]=dict(n_rows=len(r),n_galaxies=len({q['catalog_id'] for q in r}),n_exposures=len({q['filename'] for q in r}))
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
    if not all(summary.get(k,{}).get('n_rows',0)>0 for k in ['F140M','F770W']): raise ValueError('No usable overlap in at least one filter')


# Common-grid paired registration pilot. Gaussian PSF approximation only.
# Pair combinations share exposures; they are NOT independent measurements.
from scipy.ndimage import map_coordinates, gaussian_filter
STEP=.04
AX=np.arange(-50,51)*STEP
GX,GY=np.meshgrid(AX,AX)
RAD=np.hypot(GX,GY)
F140=.055
SCENARIOS=[('raw',0.,0.)]+[(f'miri{m:.2f}_target{t:.2f}',m,t)
                           for m in (.23,.26,.30) for t in (m,.40)]
def common_stamp(a,e,q,x,y,jac):
    d=np.linalg.solve(jac,np.stack([GX.ravel(),GY.ravel()])).reshape(2,*GX.shape)
    coords=[y+d[1],x+d[0]]
    # Slice detector before constructing floating masks; never interpolate bad data.
    x0=max(0,int(np.floor(np.min(coords[1])))-2);x1=min(a.shape[1],int(np.ceil(np.max(coords[1])))+3)
    y0=max(0,int(np.floor(np.min(coords[0])))-2);y1=min(a.shape[0],int(np.ceil(np.max(coords[0])))+3)
    if x1<=x0 or y1<=y0: raise ValueError('off detector')
    z=np.asarray(a[y0:y1,x0:x1],float);er=np.asarray(e[y0:y1,x0:x1],float);dq=q[y0:y1,x0:x1]
    good=np.isfinite(z)&np.isfinite(er)&(er>0)&((dq&3)==0)
    xy=[coords[0]-y0,coords[1]-x0]
    def interp(v): return map_coordinates(v,xy,order=1,mode='constant',cval=0.,prefilter=False)
    support=interp(good.astype(float))
    val=interp(np.where(good,z,0.))
    return np.divide(val,support,out=np.zeros_like(val),where=support>0),support>=.80
def blur_stamp(z,good,source,target):
    if target<=source+1e-10: return z.copy(),good.copy()
    sigma=np.sqrt(target**2-source**2)/(2.354820045*STEP)
    den=gaussian_filter(good.astype(float),sigma,mode='constant',cval=0.)
    num=gaussian_filter(np.where(good,z,0.),sigma,mode='constant',cval=0.)
    return np.divide(num,den,out=np.zeros_like(num),where=den>0),den>=.80
def plane_centroid(z,mask,r):
    ap=RAD<=r;ann=(RAD>=.95)&(RAD<=1.45)
    if mask[ap].mean()<.85 or mask[ann].mean()<.65: raise ValueError('common mask support')
    X=np.stack([np.ones_like(GX),GX,GY],axis=-1)
    M=X[ann&mask];v=z[ann&mask]
    for _ in range(5):
        b=np.linalg.lstsq(M,v,rcond=None)[0];res=v-M@b
        s=max(1e-12,1.4826*np.median(abs(res-np.median(res))))
        keep=abs(res-np.median(res))<3*s
        if keep.sum()<20 or keep.all(): break
        M=M[keep];v=v[keep]
    w=z-X@b;good=ap&mask;fl=w[good].sum()
    if not np.isfinite(fl) or fl<=0: raise ValueError('nonpositive flux')
    cen=np.array([(w[good]*GX[good]).sum(),(w[good]*GY[good]).sum()])/fl
    if np.linalg.norm(cen)>r: raise ValueError('centroid outside aperture')
    # Proxy only; interpolation/convolution correlates samples.
    sky=w[ann&mask];noise=max(1e-12,1.4826*np.median(abs(sky-np.median(sky))))
    return cen,fl/(noise*good.sum()),w
def synthetic_check():
    good=np.ones_like(GX,dtype=bool)
    a=np.exp(-((GX-.08)**2+(GY+.04)**2)/(2*.06**2))
    z,m=blur_stamp(a,good,F140,.26)
    cen,_,_=plane_centroid(z,m,.7)
    err=np.linalg.norm(cen-[.08,-.04])*1000
    print('PSF synthetic centroid error mas:',err,flush=True)
    if err>1: raise RuntimeError('synthetic centroid check failed')
def psf_pilot():
    p=argparse.ArgumentParser();p.add_argument('--refs',required=True);p.add_argument('--cal',required=True);p.add_argument('--out',required=True);args=p.parse_args()
    synthetic_check()
    out=Path(args.out)
    refs=[r for r in csv.DictReader(open(args.refs)) if r['catalog_id'] in IDS and r['rank']=='1' and float(r['separation_arcsec'])<.3]
    stamps={}; errors=[]
    for fn in sorted(glob.glob(str(Path(args.cal)/'*_cal.fits'))):
        with fits.open(fn,memmap=False) as h,datamodels.open(fn) as m:
            a,e,q=h['SCI'].data,h['ERR'].data,h['DQ'].data;w=m.meta.wcs;filt=h[0].header['FILTER']
            if filt not in ('F140M','F770W'): continue
            t=exposure_mjd(h[0].header,h['SCI'].header)
            for ref in refs:
                gid=ref['catalog_id']
                try:
                    pos=SkyCoord(float(ref['catalog_ra_deg']),float(ref['catalog_dec_deg']),unit='deg')
                    x,y=map(float,w.invert(pos.ra.deg,pos.dec.deg))
                    if not np.isfinite(x+y): raise ValueError('nonfinite WCS')
                    jac=np.zeros((2,2))
                    for j in range(2):
                        d=np.zeros(2);d[j]=.05
                        s1=SkyCoord(*map(float,w(x+d[0],y+d[1])),unit='deg');s0=SkyCoord(*map(float,w(x-d[0],y-d[1])),unit='deg')
                        dx,dy=s0.spherical_offsets_to(s1);jac[:,j]=[dx.arcsec/.1,dy.arcsec/.1]
                    z,g=common_stamp(a,e,q,x,y,jac)
                    stamps.setdefault((gid,filt),[]).append((Path(fn).name,t,z,g))
                except Exception as ex: errors.append(dict(galaxy=gid,filename=Path(fn).name,error=str(ex)))
    rows=[]; plots=set()
    for gid in sorted(IDS):
        for nf,nt,nz,ng in stamps.get((gid,'F140M'),[]):
            for mf,mt,mz,mg in stamps.get((gid,'F770W'),[]):
                # Use identical masks for raw and each matched variant of the pair.
                variants=[]
                for label,m,t in SCENARIOS:
                    if label=='raw': n,gn,k,gk=nz,ng,mz,mg
                    else:
                        n,gn=blur_stamp(nz,ng,F140,t);k,gk=blur_stamp(mz,mg,m,t)
                    variants.append((label,n,gn,k,gk))
                rawmask=ng&mg
                mask=np.logical_and.reduce([g for _,_,gn,_,gk in variants for g in (gn,gk)])
                for label,n,gn,k,gk in variants:
                    for r in (.35,.5,.7):
                        rec=dict(galaxy=gid,nircam=nf,miri=mf,nircam_mjd=nt,miri_mjd=mt,scenario=label,radius_arcsec=r,status='error',stationary_verified=False,common_ap_fraction=float(mask[RAD<=r].mean()),common_ann_fraction=float(mask[(RAD>=.95)&(RAD<=1.45)].mean()))
                        try:
                            cn,sn,wn=plane_centroid(n,mask,r);cm,sm,wm=plane_centroid(k,mask,r)
                            rec.update(status='diagnostic',east_mas=float((cm-cn)[0]*1000),north_mas=float((cm-cn)[1]*1000),nir_east_mas=float(cn[0]*1000),nir_north_mas=float(cn[1]*1000),miri_east_mas=float(cm[0]*1000),miri_north_mas=float(cm[1]*1000),nir_mean_to_sky_mad=sn,miri_mean_to_sky_mad=sm)
                            if r==.5 and label=='miri0.26_target0.26' and gid not in plots:
                                fig,axs=plt.subplots(1,2,figsize=(8,4))
                                for ax,v,c,title in zip(axs,[wn,wm],[cn,cm],['F140M blurred','F770W native']):
                                    v=np.where(mask,v,np.nan);lo,hi=np.nanpercentile(v,[5,99]);ax.imshow(v,origin='lower',extent=[AX[0],AX[-1],AX[0],AX[-1]],vmin=lo,vmax=hi,cmap='gray');ax.plot(*c,'r+');ax.add_patch(plt.Circle((0,0),r,fill=False,color='cyan'));ax.set_title(title)
                                fig.suptitle(gid+' first usable pair; Gaussian pilot');fig.tight_layout();fig.savefig(out/f'psf_pair_{gid}.png',dpi=120);plt.close(fig);plots.add(gid)
                        except Exception as ex: rec['error']=str(ex)
                        rows.append(rec)
    if not rows: raise RuntimeError('No overlapping pairs')
    with (out/'psf_matched_pair_offsets.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=sorted(set().union(*(r.keys() for r in rows))));writer.writeheader();writer.writerows(rows)
    (out/'psf_stamp_errors.json').write_text(json.dumps(errors,indent=2))
    summary=[]
    # Compare identical exposure pairs across all scenarios per galaxy and radius.
    # Equal weight per galaxy; pair medians are descriptive, not independent errors.
    for r in (.35,.5,.7):
        rr=[q for q in rows if q['radius_arcsec']==r and q['status']=='diagnostic']
        by={}
        for q in rr: by.setdefault((q['galaxy'],q['nircam'],q['miri']),{})[q['scenario']]=q
        pairs={k:v for k,v in by.items() if len(v)==len(SCENARIOS)}
        for label,_,_ in SCENARIOS:
            gal={}
            for (gid,nf,mf),v in pairs.items(): gal.setdefault(gid,[]).append([v[label]['east_mas'],v[label]['north_mas']])
            per={k:np.median(v,axis=0).tolist() for k,v in gal.items()}
            for selection in ('all','exclude_638656_sensitivity'):
                d={k:v for k,v in per.items() if selection=='all' or k!='638656'}
                if len(d)<2: continue
                a=np.array(list(d.values()));loo=np.array([np.median(np.delete(a,i,axis=0),axis=0) for i in range(len(a))])
                summary.append(dict(radius=r,scenario=label,selection=selection,n_galaxies=len(d),n_pairs=sum(len(gal[k]) for k in d),median_mas=np.median(a,axis=0).tolist(),galaxy_sd_mas=np.std(a,axis=0,ddof=1).tolist(),loo_min_mas=loo.min(axis=0).tolist(),loo_max_mas=loo.max(axis=0).tolist(),per_galaxy_mas=d))
    (out/'psf_matched_summary.json').write_text(json.dumps(summary,indent=2))
    if not summary: raise RuntimeError('No common complete scenario sample')
    for s in summary:
        if s['selection']=='exclude_638656_sensitivity' and s['scenario'] in ('raw','miri0.26_target0.26','miri0.26_target0.40'):
            print('PSF_RESULT '+json.dumps(s),flush=True)
    print('PSF pairs/rows/summary:',len({(q['galaxy'],q['nircam'],q['miri']) for q in rows}),len(rows),len(summary),flush=True)

if __name__=='__main__':
    main()
    psf_pilot()
