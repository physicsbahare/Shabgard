"""STPSF detector-sampled pilot on actual CAL cutouts; not final source association."""
import argparse, csv, glob, json
from pathlib import Path
import numpy as np
from scipy.ndimage import shift, maximum_filter, gaussian_filter
from scipy.optimize import least_squares, minimize
import stpsf

def fit(a,e,good,psf,x,y,neighbors=(),region_center=None,return_cost=False,positive=False):
    yy,xx=np.indices(a.shape)
    cx,cy=(x,y) if region_center is None else region_center
    region=good & (np.hypot(xx-cx,yy-cy)<15)
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
    lo=np.full(len(b),-np.inf)
    if positive: lo[0]=0; b[0]=max(b[0],1e-12)
    res=least_squares(lambda q:X@q-z,b,bounds=(lo,np.full(len(b),np.inf)),loss="soft_l1",f_scale=2)
    cov=np.linalg.pinv(X.T@X)
    result=(float(res.x[0]),float(np.sqrt(cov[0,0])),template(x,y))
    if return_cost: return (*result,float(2*res.cost))
    return result

def localize(a,e,good,psf,x,y,neighbors,search_radius=2.5):
    """Profile robust objective over centroid; fixed pixel mask prevents edge bias.
    Search radius is a diagnostic assumption, not a calibrated astrometric prior.
    """
    cache={}
    def objective(q):
        dx,dy=map(float,q)
        if np.hypot(dx,dy)>search_radius+1e-9:
            return 1e30
        key=(round(dx,7),round(dy,7))
        if key not in cache:
            cache[key]=fit(a,e,good,psf,x+dx,y+dy,neighbors,
                           region_center=(x,y),return_cost=True,positive=True)
        return cache[key][3]
    grid=[(dx,dy) for dx in np.linspace(-2,2,5)
          for dy in np.linspace(-2,2,5) if np.hypot(dx,dy)<=search_radius]
    seed=min(grid,key=objective)
    # SLSQP with a finite-difference step large enough for the pixel spline.
    opt=minimize(objective,seed,method="SLSQP",
                 bounds=[(-search_radius,search_radius)]*2,
                 constraints=[{"type":"ineq","fun":lambda q:search_radius**2-np.dot(q,q)}],
                 options={"eps":0.002,"maxiter":35,"ftol":1e-5})
    q=opt.x if opt.success and objective(opt.x)<=objective(seed) else np.array(seed)
    dx,dy=map(float,q)
    f,s,t,cost=fit(a,e,good,psf,x+dx,y+dy,neighbors,
                   region_center=(x,y),return_cost=True,positive=True)
    forced=objective((0.,0.))
    return dict(dx=dx,dy=dy,flux=f,sigma=s,cost=cost,
                improvement=forced-cost,optimizer_success=bool(opt.success),
                near_boundary=bool(np.hypot(dx,dy)>0.9*search_radius))

def main():
    p=argparse.ArgumentParser();p.add_argument("--root",default="miri_results");p.add_argument("--extended-controls",action="store_true");args=p.parse_args()
    root=Path(args.root)
    rows=list(csv.DictReader((root/"miri_cal_detector_diagnostic.csv").open()))
    results=[]; trials=[]
    for r in rows:
        path=root/"cutouts"/(r["filename"].replace(".fits","_local.npz"))
        if not path.exists(): continue
        d=np.load(path);a=d["sci"].astype(float);e=d["err"].astype(float)
        good=np.isfinite(a)&np.isfinite(e)&(e>0)&((d["dq"]&1)==0)
        x,y=float(d["x_target"]),float(d["y_target"])
        out={"filename":r["filename"],"filter":r["filter"],"status":"error"}
        try:
            if r.get("sci_bunit","").strip().lower().replace(" ","") not in ("mjy/sr","mjy/sr-1"):
                raise ValueError("SCI BUNIT is not validated MJy/sr")
            area=float(r["pixel_area_sr"])
            if not np.isfinite(area) or area<=0: raise ValueError("Invalid PIXAR_SR")
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
            loc=localize(a,e,good,psf,x,y,neighbors)
            null=[]; rec=[]; free_null=[]; mismatch_rec=[]; mismatch_pos=[]
            broad=gaussian_filter(psf,0.45); broad*=psf.sum()/broad.sum()
            k=area*1e12
            # Inject multiple subpixel phases and mismatch widths; preserve PSF integral.
            widths=(0.,0.25,0.45,0.65) if args.extended_controls else (0.45,)
            models={w: gaussian_filter(psf,w) if w>0 else psf.copy() for w in widths}
            for w in widths: models[w]*=psf.sum()/models[w].sum()
            for theta in np.linspace(0,2*np.pi,12,endpoint=False):
                u,v=x+18*np.cos(theta),y+18*np.sin(theta)
                try:
                    n,ns,nt=fit(a,e,good,psf,u,v,neighbors)
                    inj,_,_=fit(a+max(abs(f),5*s)*nt,e,good,psf,u,v,neighbors)
                    null.append(n);rec.append((inj-n)/max(abs(f),5*s))
                    # Same centroid search for target and controls; exclude templates coincident with trial source.
                    nn=[q for q in neighbors if np.hypot(q[0]-u,q[1]-v)>4]
                    nn.append((x+loc['dx'],y+loc['dy']))
                    z=localize(a,e,good,psf,u,v,nn)
                    free_null.append(z['flux'])
                    amp=max(loc['flux'],5*loc['sigma'])
                    _,_,template=fit(a,e,good,broad,u+0.37,v-0.29,nn)
                    zz=localize(a+amp*template,e,good,psf,u,v,nn)
                    mismatch_rec.append((zz['flux']-z['flux'])/amp)
                    mismatch_pos.append(np.hypot(zz['dx']-0.37,zz['dy']+0.29))
                    if args.extended_controls:
                        # Fit both null and injections using the identical localize routine.
                        # The same physical blank is reused across widths: not independent trials.
                        phase_x=0.45*np.cos(2*theta+0.3);phase_y=0.45*np.sin(2*theta+0.3)
                        for w in widths:
                            try:
                                _,_,tmpl=fit(a,e,good,models[w],u+phase_x,v+phase_y,nn)
                                z2=localize(a+amp*tmpl,e,good,psf,u,v,nn)
                                trials.append(dict(filename=r['filename'],filter=r['filter'],theta_rad=float(theta),
                                    blur_sigma_pix=w,phase_x_pix=float(phase_x),phase_y_pix=float(phase_y),
                                    injected_ujy=amp*k,null_free_ujy=z['flux']*k,null_formal_ujy=z['sigma']*k,
                                    null_dx_pix=z['dx'],null_dy_pix=z['dy'],null_optimizer_success=z['optimizer_success'],
                                    null_near_boundary=z['near_boundary'],injected_fit_ujy=z2['flux']*k,
                                    recovery_fraction=(z2['flux']-z['flux'])/amp,
                                    centroid_dx_error_pix=z2['dx']-phase_x,centroid_dy_error_pix=z2['dy']-phase_y,
                                    injected_optimizer_success=z2['optimizer_success'],injected_near_boundary=z2['near_boundary'],
                                    target_free_ujy=loc['flux']*k))
                            except ValueError: pass
                except ValueError: pass
            if len(null)<6: raise ValueError("insufficient nulls")
            noise=1.4826*np.median(np.abs(np.array(null)-np.median(null)))
            k=area*1e12
            fn=np.asarray(free_null)
            if len(fn)<6: raise ValueError('insufficient matched-search controls')
            out.update(matched_null_n=len(fn),matched_null_ge_target=int(np.sum(fn>=loc['flux'])),matched_null_max_ujy=float(fn.max()*k),matched_null_median_ujy=float(np.median(fn)*k),mismatch_injection_n=len(mismatch_rec),mismatch_recovery_median=float(np.median(mismatch_rec)),mismatch_centroid_error_median_pix=float(np.median(mismatch_pos)),mismatch_blur_sigma_pix=0.45)
            out.update(free_dx_pix=loc["dx"],free_dy_pix=loc["dy"],free_flux_ujy=loc["flux"]*k,
                       free_formal_sigma_ujy=loc["sigma"]*k,free_robust_cost_improvement=loc["improvement"],
                       free_optimizer_success=loc["optimizer_success"],free_near_search_boundary=loc["near_boundary"],
                       free_search_radius_pix=2.5,free_detection_validated=False,astrometric_epoch_validated=False,
                       status="physical_psf_fixed_position_DIAGNOSTIC",flux_ujy=f*k,formal_sigma_ujy=s*k,null_sigma_ujy=noise*k,n_null=len(null),null_median_ujy=float(np.median(null))*k,injection_recovery_median=float(np.median(rec)),n_neighbors=len(neighbors),psf_sum=float(psf.sum()),sigma_proxy_ujy=max(s,noise)*k)
        except Exception as exc: out["error"]=repr(exc)
        results.append(out);print(json.dumps(out),flush=True)
    fields=sorted(set().union(*(r.keys() for r in results)))
    with (root/"miri_stpsf_diagnostic.csv").open("w") as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(results)
    if args.extended_controls:
        if not trials: raise RuntimeError("No extended injection trials saved")
        with (root/"miri_extended_controls.csv").open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(trials[0]));w.writeheader();w.writerows(trials)
        print(json.dumps({"extended_trials":len(trials),"independent_trial_count_not_asserted":True}),flush=True)
    if not any(r["status"]!="error" for r in results): raise RuntimeError("No successful STPSF fits")
if __name__=="__main__":main()
