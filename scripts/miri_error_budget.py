"""Summarize actual extended MIRI controls without treating repeated widths as independent.
Run: python scripts/miri_error_budget.py miri_stpsf_diagnostic.csv miri_extended_controls.csv out
These are diagnostic quantities, not publication-ready uncertainties.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

def mad(x):
    x=np.asarray(x,float);return float(1.4826*np.median(np.abs(x-np.median(x))))

def main():
    d=pd.read_csv(sys.argv[1]);t=pd.read_csv(sys.argv[2]);out=Path(sys.argv[3]);out.mkdir(parents=True,exist_ok=True)
    assert not t.duplicated(['filename','theta_rad','blur_sigma_pix']).any()
    assert not d.duplicated(['filename']).any()
    n=t.drop_duplicates(['filename','theta_rad'])
    q=t[t.injected_optimizer_success & ~t.injected_near_boundary & t.null_optimizer_success & ~t.null_near_boundary]
    results=[]; recovery=[]
    for f,g in d.groupby('filter'):
        null=n[n['filter']==f]; allg=t[t['filter']==f]; good=q[q['filter']==f]
        m=float(g.free_flux_ujy.median())
        results.append(dict(filter=f,n_exposures=len(g),median_flux_ujy=m,dither_sample_sd_ujy=g.free_flux_ujy.std(),
            median_formal_sigma_ujy=g.free_formal_sigma_ujy.median(),matched_null_mad_ujy=mad(null.null_free_ujy),
            n_unique_controls=len(null),n_controls_converged_interior=int((null.null_optimizer_success & ~null.null_near_boundary).sum()),
            n_controls_ge_target=int((null.null_free_ujy>=null.target_free_ujy).sum()),
            max_null_to_target_ratio=(null.null_free_ujy/null.target_free_ujy).max(),n_strict_injection_trials=len(good),
            n_all_injection_trials=len(allg),final_photometry_validated=False))
        for w,a in allg.groupby('blur_sigma_pix'):
            b=good[good.blur_sigma_pix==w]
            for name,c in [('all',a),('converged_interior',b)]:
                if len(c)==0: continue
                lo,med,hi=c.recovery_fraction.quantile([.16,.5,.84])
                recovery.append(dict(filter=f,blur_sigma_pix=w,subset=name,n_trials=len(c),p16=lo,median=med,p84=hi,
                    centroid_radial_error_median_pix=np.median(np.hypot(c.centroid_dx_error_pix,c.centroid_dy_error_pix)),
                    conditional_flux_if_median_response_applied_ujy=m/med if med>0 else np.nan))
    pd.DataFrame(results).to_csv(out/'photometry_diagnostics.csv',index=False)
    pd.DataFrame(recovery).to_csv(out/'recovery_sensitivity.csv',index=False)
    print(pd.DataFrame(results).to_string(index=False));print(pd.DataFrame(recovery).to_string(index=False))

if __name__=='__main__':main()
