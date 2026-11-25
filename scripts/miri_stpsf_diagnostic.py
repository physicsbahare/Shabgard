"""STPSF detector-sampled pilot on actual CAL cutouts; not final source association."""
import argparse, csv, glob, json
from pathlib import Path
import numpy as np
from scipy.ndimage import shift, maximum_filter
from scipy.optimize import least_squares
import stpsf

def fit(a,e,good,psf,x,y,neighbors=()):
    yy,xx=np.indices(a.shape)
    region=good & (np.hypot(xx-x,yy-y)<15)
    if region.sum()<100: raise ValueError("insufficient fit pixels")
    def template(u,v):
        out=np.zeros(a.shape)
        h,w=psf.shape
        i,j=int(round(v)),int(round(u))
        p=shift(psf,(v-i,u-j),order=3,mode="constant")
        y0,x0=i-h//2,j-w//2
        ya,yb=max(0,y0),min(a.shape[0],y0+h)
        xa,xb=max(0,x0),min(a.shape[1],x0+w)
        if ya<yb and xa<xb: out[ya:yb,xa:xb]=p[ya-y0:yb-y0,xa-x0:xb-x0]
        return out
    cols=[template(x,y),np.ones(a.shape),(xx-x)/15,(yy-y)/15]+[template(u,v) for u,v in neighbors]
    X=np.stack(cols,axis=-1)[region]/e[region,None]
    z=a[region]/e[region]
    b=np.linalg.lstsq(X,z,rcond=None)[0]
    res=least_squares(lambda q:X@q-z,b,loss="soft_l1",f_scale=2)
    cov=np.linalg.pinv(X.T@X)
    return float(res.x[0]),float(np.sqrt(cov[0,0])),template(x,y)

def main():
    p=argparse.ArgumentParser();p.add_argument("--root",default="miri_results");args=p.parse_args()
    root=Path(args.root)
    rows=list(csv.DictReader((root/"miri_cal_detector_diagnostic.csv").open()))
    results=[]
    for r in rows:
        path=root/"cutouts"/(r["filename"].replace(".fits","_local.npz"))
        if not path.exists(): continue
        d=np.load(path);a=d["sci"].astype(float);e=d["err"].astype(float)
        good=np.isfinite(a)&np.isfinite(e)&(e>0)&((d["dq"]&1)==0)
        x,y=float(d["x_target"]),float(d["y_target"])
        out={"filename":r["filename"],"filter":r["filter"],"status":"error"}
        try:
            m=stpsf.MIRI();m.filter=r["filter"]
            m.detector_position=(int(round(float(r["x_photometry"]))),int(round(float(r["y_photometry"]))))
            hd=m.calc_psf(fov_pixels=41,oversample=4)
            psf=np.asarray(hd["DET_DIST"].data,dtype=float)
            # Preserve STPSF normalization; finite FOV missing wings remain systematic.
            yy,xx=np.indices(a.shape)
            bg=np.nanmedian(a[good]); sn=np.where(good,(a-bg)/e,-np.inf)
            peaks=(sn==maximum_filter(sn,size=5))&(sn>5)&(np.hypot(xx-x,yy-y)>4)&(np.hypot(xx-x,yy-y)<20)
            neighbors=[(float(u),float(v)) for v,u in np.argwhere(peaks)]
            neighbors=sorted(neighbors,key=lambda q: -sn[int(q[1]),int(q[0])])[:8]
            f,s,t=fit(a,e,good,psf,x,y,neighbors)
            null=[]; rec=[]
            for theta in np.linspace(0,2*np.pi,12,endpoint=False):
                u,v=x+18*np.cos(theta),y+18*np.sin(theta)
                try:
                    n,ns,nt=fit(a,e,good,psf,u,v,neighbors)
                    inj,_,_=fit(a+max(abs(f),5*s)*nt,e,good,psf,u,v,neighbors)
                    null.append(n);rec.append((inj-n)/max(abs(f),5*s))
                except ValueError: pass
            if len(null)<6: raise ValueError("insufficient nulls")
            noise=1.4826*np.median(np.abs(np.array(null)-np.median(null)))
            k=float(r["pixel_area_sr"])*1e12
            out.update(status="physical_psf_fixed_position_DIAGNOSTIC",flux_ujy=f*k,formal_sigma_ujy=s*k,null_sigma_ujy=noise*k,n_null=len(null),null_median_ujy=float(np.median(null))*k,injection_recovery_median=float(np.median(rec)),n_neighbors=len(neighbors),psf_sum=float(psf.sum()),sigma_proxy_ujy=max(s,noise)*k)
        except Exception as exc: out["error"]=repr(exc)
        results.append(out);print(json.dumps(out),flush=True)
    fields=sorted(set().union(*(r.keys() for r in results)))
    with (root/"miri_stpsf_diagnostic.csv").open("w") as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(results)
    if not any(r["status"]!="error" for r in results): raise RuntimeError("No successful STPSF fits")
if __name__=="__main__":main()
