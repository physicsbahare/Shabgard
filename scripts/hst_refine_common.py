"""Refine HST target likelihood on original FLTs with fixed common reference subsets.
Reuse archived reference measurements; compare subsets and PSF widths, not independent epochs.
"""
import csv,json,glob
from collections import defaultdict
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import hst_wfc3_ir_epoch as h
from pathlib import Path
def read(n):
    p=glob.glob('pilot/**/'+n,recursive=True)
    assert len(p)==1,p
    return list(csv.DictReader(open(p[0])))
def main():
    rows=read('exposure_audit.csv'); refs=defaultdict(dict)
    for q in read('references.csv'):
        refs[q['filename']][str(q['ref_id'])]=[float(q['east_mas']),float(q['north_mas'])]
    out=Path('hst_refined');out.mkdir(exist_ok=True)
    groups=defaultdict(list); measurements=[]
    for q in rows:
        if q['status']!='pixel_likelihood':continue
        n=q['filename'];date=q['visit'];f=q['filter'];t=float(q['mid_mjd'])
        p=h.download(n,'mast:HST/product/'+n)
        with fits.open(p,memmap=False) as hd:
            sh=hd['SCI',1];a=sh.data;e=hd['ERR',1].data;d=hd['DQ',1].data;w=WCS(sh.header,hd)
            pos=h.pred(t)
            subsets=[('common4',['6','10','11','16']),('common3',['6','10','11'])] if date!='2010-11-03' else [('unregistered',[])]
            for label,ids in subsets:
                if not all(k in refs[n] for k in ids):raise ValueError((n,ids,list(refs[n])))
                c=np.median([refs[n][k] for k in ids],axis=0) if ids else np.zeros(2)
                center=pos.spherical_offsets_by(c[0]*u.mas,c[1]*u.mas)
                # Broad enough to cover previous maxima, 5 mas grid.
                xy=np.arange(-.20,.201,.005);xx,yy=np.meshgrid(xy,xy);grid=np.c_[xx.ravel(),yy.ravel()]
                for width in [.85,1.,1.15]:
                    sig={'F125W':.55,'F140W':.60,'F160W':.65}[f]*width
                    ff,vv=h.scan(a,e,d,w,center,grid,sig);sn=h.combined([(ff,vv)])
                    i=np.nanargmax(sn);key=(date,f,label,width)
                    groups[key].append((ff,vv))
                    measurements.append(dict(filename=n,visit=date,filter=f,subset=label,width=width,mid_mjd=t,dx_mas=grid[i,0]*1000,dy_mas=grid[i,1]*1000,snr=float(sn[i]),boundary=bool(np.max(abs(grid[i]))>.195),translation_east_mas=c[0],translation_north_mas=c[1]))
        print('REFINED_EXPOSURE',n,flush=True)
    summary=[]
    for key,g in groups.items():
        sn=h.combined(g);i=np.nanargmax(sn);date,f,label,width=key
        # Profile likelihood curvature for positive amplitude, not full systematic covariance.
        delta=grid-grid[i];sel=(np.linalg.norm(delta,axis=1)<=.025)&np.isfinite(sn)
        x,y=delta[sel].T*1000
        X=np.c_[np.ones(len(x)),x,y,x*x/2,x*y,y*y/2]
        b=np.linalg.lstsq(X,-.5*sn[sel]**2,rcond=None)[0]
        H=np.array([[b[3],b[4]],[b[4],b[5]]])
        cov=np.linalg.inv(H) if np.min(np.linalg.eigvalsh(H))>0 else None
        rec=dict(visit=date,filter=f,subset=label,width=width,n=len(g),dx_mas=grid[i,0]*1000,dy_mas=grid[i,1]*1000,snr=float(sn[i]),boundary=bool(np.max(abs(grid[i]))>.195),formal_cov_mas2=cov.tolist() if cov is not None else None)
        summary.append(rec);print('REFINED_VISIT',json.dumps(rec),flush=True)
        np.savez_compressed(out/f'{date}_{f}_{label}_{width}.npz',grid=grid,snr=sn)
    for name,data in [('exposures',measurements),('summary',summary)]:
        (out/(name+'.json')).write_text(json.dumps(data,indent=2))
    assert summary
if __name__=='__main__':main()
