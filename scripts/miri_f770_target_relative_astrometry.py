"""Map actual free-centroid STPSF offsets to sky with original F770W CAL gWCS.
Relative frame only: no external galaxy tie or parallax inference.
"""
import argparse,csv,json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from jwst import datamodels
from audit_calibrated_detector_local import exposure_mjd
from scipy.ndimage import map_coordinates
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def load(path):
    return list(csv.DictReader(open(path)))

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--cal',required=True);p.add_argument('--diagnostic',required=True)
    p.add_argument('--psf',required=True);p.add_argument('--network',required=True)
    p.add_argument('--out',required=True);args=p.parse_args()
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    diag={r['filename']:r for r in load(args.diagnostic) if r.get('filter')=='F770W'}
    psf={r['filename']:r for r in load(args.psf) if r.get('filter')=='F770W'}
    netroot=Path(args.network)
    net=next(r for r in json.loads((netroot/'network.json').read_text()) if r['radius_mas']==250 and r['snr_cut']==0)
    cv=np.load(netroot/'relative_translation_covariance.npz',allow_pickle=False)
    names=list(cv['filenames']); cov=cv['covariance_mas2']
    assert cov.shape==(2*len(names),2*len(names))
    offsets={r['filename']:r for r in net['offsets']}
    if set(names)!=set(offsets): raise ValueError('Covariance/offset file mismatch')
    center=SkyCoord(150.0634341,2.3552216,unit='deg')
    rows=[];stamps=[];vectors=[]
    for filename in names:
        r=diag[filename];q=psf[filename];o=offsets[filename]
        if q['free_optimizer_success'].lower()!='true' or q['free_near_search_boundary'].lower()=='true':
            raise ValueError('Unvalidated centroid search '+filename)
        x=float(r['x_photometry'])+float(q['free_dx_pix'])
        y=float(r['y_photometry'])+float(q['free_dy_pix'])
        with fits.open(Path(args.cal)/filename,memmap=False) as hd,datamodels.open(str(Path(args.cal)/filename)) as model:
            w=model.meta.wcs
            ra,dec=map(float,w(x,y))
            if not np.isfinite(ra+dec): raise ValueError('Nonfinite gWCS coordinates')
            sky=SkyCoord(ra,dec,unit='deg')
            east,north=center.spherical_offsets_to(sky)
            raw=np.array([east.mas,north.mas]); corr=raw-np.array([o['ra_mas'],o['dec_mas']])
            csky=center.spherical_offsets_by(corr[0]*u.mas,corr[1]*u.mas)
            jac=np.zeros((2,2))
            for j in range(2):
                d=np.zeros(2);d[j]=0.01
                plus=SkyCoord(*map(float,w(x+d[0],y+d[1])),unit='deg')
                minus=SkyCoord(*map(float,w(x-d[0],y-d[1])),unit='deg')
                a,b=minus.spherical_offsets_to(plus);jac[:,j]=[a.mas/0.02,b.mas/0.02]
            if not np.all(np.isfinite(jac)): raise ValueError('Bad gWCS Jacobian')
            phead=hd[0].header;sh=hd['SCI'].header
            t=exposure_mjd(phead,sh)
            dq=np.asarray(hd['DQ'].data);a=np.asarray(hd['SCI'].data,float)
            if not (10<x<a.shape[1]-10 and 10<y<a.shape[0]-10): raise ValueError('Centroid outside CAL safe margin')
            gy,gx=np.mgrid[-9:10,-9:10]
            stamp=map_coordinates(np.where((dq&1)==0,a,np.nan),[y+gy,x+gx],order=1,mode='constant',cval=np.nan)
            stamps.append(stamp);vectors.append([raw,corr])
            rec=dict(filename=filename,mid_mjd=t,x_detector=x,y_detector=y,ra_raw_deg=ra,dec_raw_deg=dec,
                relative_ra_deg=csky.ra.deg,relative_dec_deg=csky.dec.deg,raw_east_mas=raw[0],raw_north_mas=raw[1],
                relative_east_mas=corr[0],relative_north_mas=corr[1],translation_east_mas=o['ra_mas'],translation_north_mas=o['dec_mas'],
                jac_east_x=jac[0,0],jac_east_y=jac[0,1],jac_north_x=jac[1,0],jac_north_y=jac[1,1],
                free_flux_ujy=float(q['free_flux_ujy']),absolute_frame_validated=False,target_centroid_covariance_validated=False,
                mjd_avg_header=phead.get('MJD-AVG',sh.get('MJD-AVG','')),
                expstart=phead.get('EXPSTART',''),expend=phead.get('EXPEND',''))
            for key in ('JWST_X','JWST_Y','JWST_Z','JWST_DX','JWST_DY','JWST_DZ','EPH_TIME','TIMESYS','TREFPOS'):
                rec[key]=phead.get(key,sh.get(key,''))
            rows.append(rec)
    with (out/'target_relative_positions.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    v=np.array(vectors);n=len(rows)
    A=np.tile(np.eye(2)/n,(1,n))
    prop=A@cov@A.T
    # Gauge-only, illustrative equal weight; never a final target position error.
    summary=dict(n_exposures=n,raw_mean_mas=v[:,0].mean(axis=0).tolist(),relative_mean_mas=v[:,1].mean(axis=0).tolist(),
        raw_sample_sd_mas=v[:,0].std(axis=0,ddof=1).tolist(),relative_sample_sd_mas=v[:,1].std(axis=0,ddof=1).tolist(),
        translation_only_mean_covariance_mas2=prop.tolist(),relative_frame_reference=net['reference'],
        absolute_frame_validated=False,parallax_ready=False,
        note='Target centroid noise, PSF mismatch, external frame and chromatic offsets not included; covariance only describes relative translations.')
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)
    fig,axs=plt.subplots(2,4,figsize=(12,6))
    for ax,s,r in zip(axs.flat,stamps,rows):
        lo,hi=np.nanpercentile(s,[10,99]);ax.imshow(s,origin='lower',vmin=lo,vmax=hi,cmap='gray')
        ax.plot(9,9,'r+',ms=10);ax.set_title(r['filename'].split('_mirimage')[0][-11:]);ax.set_axis_off()
    fig.tight_layout();fig.savefig(out/'target_cal_centroid_stamps.png',dpi=160);plt.close(fig)
    np.savez(out/'target_relative_translation_covariance.npz',covariance_mas2=cov,filenames=np.array(names),axis_order=cv['axis_order'])

if __name__=='__main__':main()
