"""Common-reference jackknife on actual HST subpixel and registration artifacts.
All exposures share the same deleted reference; retain cross-visit covariance.
Actual target ePSF centroids from run 46121176137 logs; reference registration remains diagnostic."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
def locate(root,name):
    p=list(Path(root).rglob(name))
    assert len(p)==1,(root,name,p)
    return p[0]
a=pd.read_csv(locate('hst','exposure_audit.csv'))
r=pd.read_csv(locate('hst','references.csv'))
s=pd.DataFrame(json.loads("[{\"filename\":\"ibhm55e7q_flt.fits\",\"east\":124.6336361370673,\"north\":38.467097958184766,\"valid\":true},{\"filename\":\"ibhm55eeq_flt.fits\",\"east\":107.32450155632898,\"north\":59.76089866392964,\"valid\":true},{\"filename\":\"ibhm55elq_flt.fits\",\"east\":83.19117960848797,\"north\":44.37597868701511,\"valid\":true},{\"filename\":\"ibhm55esq_flt.fits\",\"east\":107.93448685669051,\"north\":44.52499397118092,\"valid\":true},{\"filename\":\"iboa30rlq_flt.fits\",\"east\":0.11121577738581594,\"north\":18.6803980546331,\"valid\":true},{\"filename\":\"iboa30rnq_flt.fits\",\"east\":-1.3420953511597362,\"north\":24.25776345494768,\"valid\":true},{\"filename\":\"iboa30rrq_flt.fits\",\"east\":-6.431066970989353,\"north\":17.51615047232101,\"valid\":true},{\"filename\":\"iboa30ruq_flt.fits\",\"east\":10.66870493253689,\"north\":15.125025585259863,\"valid\":true},{\"filename\":\"iboa74bmq_flt.fits\",\"east\":0.7851686721175417,\"north\":19.496233447885547,\"valid\":true},{\"filename\":\"iboa74bpq_flt.fits\",\"east\":-15.008169296197728,\"north\":10.50670705655768,\"valid\":true},{\"filename\":\"iboa74bsq_flt.fits\",\"east\":-8.583153785944734,\"north\":27.901884132787934,\"valid\":true},{\"filename\":\"iboa74bwq_flt.fits\",\"east\":-12.518646860377675,\"north\":30.688353938612295,\"valid\":true}]"))
s=s[s.valid==True]
a=a[a.registered==True].merge(s,on='filename',validate='one_to_one')
a=a.sort_values(['visit','filter','filename'])
assert len(a)==8,len(a)
ids=sorted(set.intersection(*(set(r.loc[r.filename==f,'ref_id']) for f in a.filename)))
assert len(ids)>=3,ids
z=np.stack([r[r.filename==f].set_index('ref_id').loc[ids,['east_mas','north_mas']].to_numpy(float) for f in a.filename])
base=a[['east','north']].to_numpy()+a[['translation_east_mas','translation_north_mas']].to_numpy()
visits=sorted(a.visit.unique())
def calc(k):
    xy=base-np.median(z[:,k,:],axis=1)
    # Equal filter weighting per visit, independent of exposure SNR.
    filt=np.array([np.mean(xy[(a.visit==v)&(a['filter']==f)],axis=0) for v in visits for f in ['F125W','F160W']])
    vv=filt.reshape(len(visits),2,2).mean(axis=1)
    return xy,filt,vv
xy,filt,v=calc(np.arange(len(ids)))
loo=np.array([calc(np.delete(np.arange(len(ids)),i))[2].ravel() for i in range(len(ids))])
d=loo-loo.mean(axis=0)
C=(len(ids)-1)/len(ids)*d.T@d
res=dict(ids=[int(i) for i in ids],visits=visits,visit_offsets_mas=v.tolist(),filter_differences_mas=(filt[::2]-filt[1::2]).tolist(),leave_one_reference_out_mas=loo.tolist(),shared_reference_covariance_mas2=C.tolist(),marginal_sigma_mas=np.sqrt(np.diag(C)).tolist(),rank=int(np.linalg.matrix_rank(C)),within_visit_exposure_scatter_mas=[np.std(xy[a.visit==t],axis=0,ddof=1).tolist() for t in visits])
print('REFERENCE_JACKKNIFE',json.dumps(res))
Path('result').mkdir(exist_ok=True)
Path('result/reference_jackknife.json').write_text(json.dumps(res,indent=2))
a[['filename','visit','filter']].assign(east=xy[:,0],north=xy[:,1]).to_csv('result/reregistered_exposures.csv',index=False)
