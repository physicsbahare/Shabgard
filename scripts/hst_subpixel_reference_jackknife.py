"""Common-reference jackknife on actual HST subpixel and registration artifacts.
All exposures share the same deleted reference; retain cross-visit covariance.
Diagnostic only, not absolute calibration."""
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
s=pd.read_csv(locate('subpixel','subpixel.csv'))
s=s[s.valid==True].copy()
s['filename']=s.file.str.replace('_likelihood.npz','',regex=False)
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
