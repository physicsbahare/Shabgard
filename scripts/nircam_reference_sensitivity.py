"""Reproduce correlated reference resampling on actual archived JWST measurements.
Input ZIP names and artifact IDs are in docs/reference_sensitivity_20261127.md.
The bootstrap is conditional on successful registration, not a calibrated systematic floor.
"""
import numpy as np,pandas as pd,json,zipfile,sys
from pathlib import Path
root=Path(sys.argv[1] if len(sys.argv)>1 else '.')
def read(z,n):
 with zipfile.ZipFile(root/z) as f:return pd.read_csv(f.open(n))
def js(z,n):
 with zipfile.ZipFile(root/z) as f:return json.load(f.open(n))
b=read('baseline_refaudit.zip','reference_centroids.csv');r=read('baseline_refaudit.zip','registration_pilot.csv');t=read('baseline_refaudit.zip','target_centroids.csv');f=read('fourth_refaudit.zip','references.csv');ft=read('fourth_refaudit.zip','targets.csv');g=js('gaia_checked.zip','matches.json');j=js('joint_exact.zip','joint_pilot_result.json');fac=read('joint_exact.zip','hst_factors.csv')
ids=sorted(b.groupby('ref_id').program.nunique().loc[lambda x:x==3].index);b=b[b.ref_id.isin(ids)].copy()
def trans(v,minimum):
 if len(v)<minimum:return None
 c=np.median(v,axis=0);mad=np.maximum(1,1.4826*np.median(abs(v-c),axis=0));v=v[np.all(abs(v-c)<=4*mad,axis=1)]
 return np.median(v,axis=0) if len(v)>=minimum else None
center=np.array([150.0634341,2.3552216]);scale=np.array([np.cos(np.deg2rad(center[1])),1])*3.6e6
orig=r.set_index('filename')[['east_mas','north_mas']];b=b[b.filename.isin(orig.index)].copy()
b[['old_e','old_n']]=b[['east_mas','north_mas']].to_numpy()-orig.loc[b.filename].to_numpy()
old=b.groupby(['ref_id','program'])[['old_e','old_n']].median().groupby('ref_id').median()
raw_t=t.set_index('filename');raw_f=ft.set_index('filename')
def reconstruct(sample):
 counts=pd.Series(sample).value_counts();corr={};points=[]
 for name,z in b.groupby('filename'):
  z=z[z.ref_id.isin(counts.index)];v=np.repeat(z[['east_mas','north_mas']].to_numpy(),counts.loc[z.ref_id].to_numpy(),axis=0);c=trans(v,6)
  if c is None:return None
  corr[name]=c;p=(raw_t.loc[name,['ra_deg','dec_deg']].to_numpy(float)-center)*scale-c;points.append([int(z.program.iloc[0]),*p])
 bb=b.copy();bb[['new_e','new_n']]=bb[['east_mas','north_mas']].to_numpy()-np.array([corr[n] for n in bb.filename]);base=bb.groupby(['ref_id','program'])[['new_e','new_n']].median().groupby('ref_id').median();base.columns=['old_e','old_n'];delta=base-old
 for name,z in f.groupby('filename'):
  z=z[(z.ref_id.isin(counts.index))&(z.baseline_epoch_scatter_mas<=10)];v=z[['east','north']].to_numpy()-delta.loc[z.ref_id].to_numpy();v=np.repeat(v,counts.loc[z.ref_id].to_numpy(),axis=0);c=trans(v,3)
  if c is None:return None
  p=(raw_f.loc[name,['ra','dec']].to_numpy(float)-center)*scale-c;points.append([7814,*p])
 return pd.DataFrame(points,columns=['program','e','n']).groupby('program')[['e','n']].median().loc[[1837,1727,5893,7814]].to_numpy().ravel()
basepos=reconstruct(ids);y=np.array(j['rows'])[:,3:].ravel();assert np.max(abs(basepos-y[4:]))<.1
rows=np.array(j['rows']);T=(rows[:,0]-60308.94)/365.25;A=np.zeros((12,5));A[::2,0]=1;A[1::2,1]=1;A[::2,2]=T;A[1::2,3]=T;A[:,4]=rows[:,1:3].ravel();C=np.diag(np.repeat([10,10,5,5,5,10],2)**2).astype(float);C[:2,:2]+=[[43.2292,-7.2917],[-7.2917,10.4167]];C[2:4,2:4]+=[[30.7292,25],[25,20.8333]];v=np.r_[j['hst_reference_shift'],np.zeros(8)];C+=np.outer(v,v);U=np.zeros((12,2));U[:4]=np.tile(np.eye(2),(2,1));C+=100*U@U.T
W=np.linalg.inv(C);M=np.linalg.solve(A.T@W@A,A.T@W);p=M@y
loo={}
for k in ids:
 d=reconstruct([i for i in ids if i!=k])
 if d is not None:loo[int(k)]=dict(shift_mas=(d-basepos).tolist(),delta_pi=float(M[4,4:]@(d-basepos)))
rng=np.random.default_rng(639188);boot=[]
for i in range(180):
 d=reconstruct(rng.choice(ids,len(ids),replace=True))
 if d is not None:boot.append(d-basepos)
boot=np.array(boot);cov=np.cov(boot,rowvar=False);pi_delta=boot@M[4,4:];C2=C.copy();C2[4:,4:]+=cov;W2=np.linalg.inv(C2);V2=np.linalg.inv(A.T@W2@A);p2=V2@A.T@W2@y
out=dict(nrefs=len(ids),reconstruction_delta=(basepos-y[4:]).tolist(),gaia_matches=[q for q in g if q.get('matches')],gaia_no_match=sum(not q.get('matches') for q in g),loo=loo,bootstrap_n=len(boot),bootstrap_attempts=180,bootstrap_cov_mas2=cov.tolist(),bootstrap_pi_sd=float(np.std(pi_delta,ddof=1)),bootstrap_pi_quantiles=np.quantile(pi_delta,[.025,.5,.975]).tolist(),baseline_pi=float(p[4]),baseline_sigma=float(np.sqrt(np.linalg.inv(A.T@W@A)[4,4])),illustrative_added_cov_pi=float(p2[4]),illustrative_added_cov_sigma=float(np.sqrt(V2[4,4])),hst_orbit_au_range=[float(fac.hst_earth_au.min()),float(fac.hst_earth_au.max())],max_hst_earth_factor_difference=float(np.max(np.hypot(fac.px-fac.earth_px,fac.py-fac.earth_py))))
print(json.dumps(out,indent=2));(root/'reference_sensitivity_result.json').write_text(json.dumps(out,indent=2))
