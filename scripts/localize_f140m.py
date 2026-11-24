import csv,glob,json
import numpy as np
from scipy.special import erf
from scipy.optimize import least_squares
from astropy.io import fits
from jwst import datamodels
from audit_calibrated_detector_local import exposure_mjd,sky_at

def fit(a,e,d,x,y,search=3):
    ix,iy=int(round(x)),int(round(y))
    a=a[iy-7:iy+8,ix-7:ix+8]; e=e[iy-7:iy+8,ix-7:ix+8]; d=d[iy-7:iy+8,ix-7:ix+8]
    if a.shape!=(15,15): raise ValueError("edge")
    yy,xx=np.indices(a.shape); xx=xx+ix-7-x; yy=yy+iy-7-y
    good=np.isfinite(a)&np.isfinite(e)&(e>0)&(d==0)
    if good.sum()<100: raise ValueError("masked")
    b=np.median(a[good]); floor=1.4826*np.median(abs(a[good]-b)); e=np.maximum(e,floor)
    def model(p):
        f,u,v,s,c,g,h=p
        def axis(z,t): return .5*(erf((z-t+.5)/(np.sqrt(2)*s))-erf((z-t-.5)/(np.sqrt(2)*s)))
        return f*axis(xx,u)*axis(yy,v)+c+g*xx+h*yy
    def residual(p): return ((model(p)-a)[good]/e[good])
    scale=max(floor,np.median(e[good]),1e-12)
    lo=[0,-search,-search,.45,-np.inf,-np.inf,-np.inf]; hi=[np.inf,search,search,1.4,np.inf,np.inf,np.inf]
    fits0=[least_squares(residual,[scale*10,u,v,.7,b,0,0],bounds=(lo,hi),max_nfev=150) for u,v in [(0,0),(-2,-2),(-2,2),(2,-2),(2,2)]]
    r=min(fits0,key=lambda q: q.cost)
    p=r.x; cov=np.linalg.pinv(r.jac.T@r.jac)*max(1,2*r.cost/(good.sum()-7))
    null=np.linalg.lstsq(np.stack([np.ones_like(xx),xx,yy],axis=-1)[good]/e[good,None],a[good]/e[good],rcond=None)[0]
    chi0=np.sum(((a-null[0]-null[1]*xx-null[2]*yy)[good]/e[good])**2)
    return dict(x=float(x+p[1]),y=float(y+p[2]),flux=float(p[0]),sigma=float(p[3]),formal_x=float(np.sqrt(cov[1,1])),formal_y=float(np.sqrt(cov[2,2])),delta_chi2=float(max(0,chi0-2*r.cost)),reduced_chi2=float(2*r.cost/(good.sum()-7)),boundary=bool(np.any(r.active_mask)),converged=bool(r.success))

def main():
    rows=[]
    for path in sorted(glob.glob('pilot_data/*_cal.fits')):
        with fits.open(path,memmap=False) as h, datamodels.open(path) as m:
            a,e,d=h['SCI'].data,h['ERR'].data,h['DQ'].data
            t=exposure_mjd(h[0].header,h['SCI'].header); w=m.meta.wcs
            x,y=w.invert(*sky_at(t,60308))
            q=fit(a,e,d,x,y); ra,dec=w(q['x'],q['y'])
            controls=[]
            for dx in range(-100,101,20):
                for dy in range(-100,101,20):
                    if dx==dy==0: continue
                    try: controls.append(fit(a,e,d,x+dx,y+dy)['delta_chi2'])
                    except ValueError: pass
            q.update(path=path,mjd=t,ra=float(ra),dec=float(dec),pred_x=float(x),pred_y=float(y),control_n=len(controls),control_ge=sum(z>=q['delta_chi2'] for z in controls),control_delta_chi2=controls)
            rows.append(q)
    assert len(rows)==3
    with open('pilot_results/localization.json','w') as f: json.dump(rows,f,indent=2)
    for q in rows: print(json.dumps({k:v for k,v in q.items() if k!='control_delta_chi2'}))
if __name__=='__main__': main()
