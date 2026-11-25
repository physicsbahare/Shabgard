"""F770W same-visit field-control pilot on original calibrated detector images.
Bright isolated compact *candidates*, not spectroscopically classified stars.
Measures real profile concentration relative to STPSF and repeat-dither
translation; not an absolute or cross-instrument astrometric registration.
"""
import argparse,csv,json,glob
from pathlib import Path
import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter,maximum_filter
from scipy.spatial import cKDTree
from jwst import datamodels
import stpsf
from miri_stpsf_diagnostic import localize,fit
from audit_calibrated_detector_local import exposure_mjd,sky_at

def dump(path,rows):
    if not rows: return
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=sorted(set().union(*(r.keys() for r in rows))))
        w.writeheader();w.writerows(rows)

def tangent(ra,dec,ra0,dec0):
    return np.array([(ra-ra0)*np.cos(np.deg2rad(dec0))*3600,(dec-dec0)*3600])

def measure(p):
    rows=[]
    with fits.open(p,memmap=False) as h,datamodels.open(p) as model:
        a=np.asarray(h["SCI"].data,float);e=np.asarray(h["ERR"].data,float)
        q=np.asarray(h["DQ"].data,np.uint32)
        if str(h[0].header.get("FILTER","")).upper()!="F770W": return rows
        if str(h["SCI"].header.get("BUNIT","")).strip().lower()!="mjy/sr":
            raise ValueError("unexpected SCI unit")
        w=model.meta.wcs
        t=exposure_mjd(h[0].header,h["SCI"].header)
        r0,d0=sky_at(t,60308.)
        tx,ty=map(float,w.invert(r0,d0))
        if not np.isfinite(tx+ty): raise ValueError("no target gWCS")
        good=np.isfinite(a)&np.isfinite(e)&(e>0)&((q&3)==0)
        # Filtered local excess for selecting candidates; do not use as a flux or significance.
        filled=np.where(good,a,np.nanmedian(a[good]))
        sm=gaussian_filter(filled,1.0)-gaussian_filter(filled,6.)
        sn=np.where(good,sm/e,-np.inf)
        yy,xx=np.indices(a.shape)
        rr=np.hypot(xx-tx,yy-ty)
        # Stay near the target to limit field-dependent optics/distortion.
        peaks=(sn==maximum_filter(sn, size=7))&(sn>5)&(rr>30)&(rr<180)
        peaks &= (xx>=25)&(xx<a.shape[1]-25)&(yy>=25)&(yy<a.shape[0]-25)
        xy=np.argwhere(peaks)
        xy=sorted(xy,key=lambda z:-sn[z[0],z[1]])[:50]
        m=stpsf.MIRI();m.filter="F770W"
        m.detector_position=(int(round(tx)),int(round(ty)))
        psf=np.asarray(m.calc_psf(fov_pixels=41,oversample=4)["DET_DIST"].data,float)
        # All target and field candidates use identical fitting PSF.
        for v,u in xy:
            # Require isolation from other significant local maxima.
            if any(0<np.hypot(u-uu,v-vv)<18 for vv,uu in xy): continue
            sl=np.s_[v-22:v+23,u-22:u+23]
            aa,ee,gg=a[sl],e[sl],good[sl]
            gy,gx=np.indices(aa.shape)
            if np.mean(gg[np.hypot(gx-22,gy-22)<9])<0.98: continue
            try:
                loc=localize(aa,ee,gg,psf,22.,22.,[])
                if not loc["optimizer_success"] or loc["near_boundary"] or loc["flux"]<=0: continue
                # Profile shape from data, independent of fitted amplitude.
                cx,cy=22+loc["dx"],22+loc["dy"]
                rad=np.hypot(gx-cx,gy-cy)
                ann=gg&(rad>=12)&(rad<=18)
                if ann.sum()<150: continue
                X=np.stack([np.ones(aa.shape),(gx-cx)/15,(gy-cy)/15],-1)
                z=aa[ann]; M=X[ann]
                b=np.linalg.lstsq(M,z,rcond=None)[0]
                for _ in range(3):
                    res=z-M@b; med=np.median(res)
                    s=1.4826*np.median(np.abs(res-med))
                    use=np.abs(res-med)<3*max(s,1e-12)
                    if use.sum()<50: break
                    b=np.linalg.lstsq(M[use],z[use],rcond=None)[0]
                resid=aa-X@b
                _,_,template=fit(aa,ee,gg,psf,cx,cy,[])
                core=gg&(rad<=3); outer=gg&(rad<=8)
                total=resid[outer].sum()
                if total<=0: continue
                ratio=resid[core].sum()/total
                modelratio=template[core].sum()/template[outer].sum()
                ra,dec=map(float,w(u+loc["dx"],v+loc["dy"]))
                if not np.isfinite(ra+dec): continue
                rows.append(dict(filename=Path(p).name,mid_mjd=t,x=u+loc["dx"],y=v+loc["dy"],
                    ra=ra,dec=dec,profile_core_to_8pix=ratio,stpsf_core_to_8pix=modelratio,
                    profile_ratio=ratio/modelratio,fit_snr_proxy=loc["flux"]/loc["sigma"],
                    distance_to_target_pix=float(np.hypot(u-tx,v-ty)),
                    classification="unclassified_isolated_compact_candidate"))
            except (ValueError,RuntimeError,np.linalg.LinAlgError): continue
    return rows

def main():
    ap=argparse.ArgumentParser();ap.add_argument("files",nargs="+");ap.add_argument("--out",default="field_results")
    args=ap.parse_args(); out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    paths=sorted(set(p for pat in args.files for p in glob.glob(pat)))
    allrows=[];byfile={}
    for p in paths:
        try:
            r=measure(p);byfile[Path(p).name]=r;allrows+=r
            print(json.dumps({"filename":Path(p).name,"field_candidates":len(r)}),flush=True)
        except Exception as e: print(json.dumps({"filename":p,"error":repr(e)}),flush=True)
    dump(out/"field_candidates.csv",allrows)
    if not allrows: raise SystemExit("no field candidates measured")
    ref=max(byfile,key=lambda k:len(byfile[k]))
    ra0,dec0=150.0634341,2.3552216
    def coords(rows):
        return np.array([tangent(r["ra"],r["dec"],ra0,dec0) for r in rows])
    base=coords(byfile[ref]); tree=cKDTree(base)
    matches=[];summaries=[]
    for f,r in byfile.items():
        if f==ref or not r: continue
        pos=coords(r);dist,idx=tree.query(pos)
        back=cKDTree(pos).query(base)[1]
        pairs=[(i,j) for i,j in enumerate(idx) if dist[i]<0.5 and back[j]==i]
        if len(pairs)<3: continue
        delta=np.array([(pos[i]-base[j])*1000 for i,j in pairs])
        keep=np.ones(len(delta),bool)
        for _ in range(3):
            med=np.median(delta[keep],axis=0)
            norm=np.linalg.norm(delta-med,axis=1)
            scale=1.4826*np.median(np.abs(norm[keep]-np.median(norm[keep])))
            k=norm<max(25.,np.median(norm[keep])+3*scale)
            if k.sum()<3: break
            keep=k
        med=np.median(delta[keep],axis=0); scatter=1.4826*np.median(np.abs(delta[keep]-med),axis=0)
        for n,(i,j) in enumerate(pairs):
            matches.append(dict(filename=f,reference_filename=ref,reference_index=int(j),
                delta_ra_cosdec_mas=delta[n,0],delta_dec_mas=delta[n,1],used=bool(keep[n]),
                profile_ratio=r[i]["profile_ratio"]))
        rec=dict(filename=f,reference_filename=ref,n_matches=len(pairs),n_used=int(keep.sum()),
            translation_ra_mas=med[0],translation_dec_mas=med[1],
            residual_mad_ra_mas=scatter[0],residual_mad_dec_mas=scatter[1])
        summaries.append(rec); print(json.dumps(rec),flush=True)
    dump(out/"same_visit_matches.csv",matches);dump(out/"same_visit_translation.csv",summaries)
    # Report each reference object only once, not every dither as independent.
    ids=sorted(set(r["reference_index"] for r in matches if r["used"]))
    repeated=[byfile[ref][i]["profile_ratio"] for i in ids]
    summary=dict(n_files=len(byfile),n_distinct_matched_controls=len(ids),
        n_registered_dithers=len(summaries),profile_ratio_median=float(np.median(repeated)) if repeated else None,
        profile_ratio_mad=float(1.4826*np.median(np.abs(np.array(repeated)-np.median(repeated)))) if repeated else None,
        reference_filename=ref,stationary_classification_verified=False,
        cross_epoch_registration_validated=False,photometry_publication_ready=False)
    (out/"summary.json").write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
    if not summaries: raise SystemExit("no repeated field matches; inspect candidate selection")
if __name__=="__main__": main()
