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
if __name__=='__main__':main()
