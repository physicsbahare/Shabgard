"""Six-epoch GLS pilot, unconstrained parallax; inputs are three actual Actions ZIPs.
Requires numpy, pandas; swisseph only for explicitly approximate Earth fallback.
"""
import numpy as np,pandas as pd,json,zipfile,sys
from pathlib import Path
root=Path(sys.argv[1] if len(sys.argv)>1 else '.')
def read(z,n):
 with zipfile.ZipFile(root/z) as f:return pd.read_csv(f.open(n))
def js(z,n):
 with zipfile.ZipFile(root/z) as f:return json.load(f.open(n))
r=read('nircam_baseline_fit.zip','registration_pilot.csv');f=read('nircam_fourth_fit.zip','jwst_parallax_factors.csv');q=read('nircam_fourth_fit.zip','registration.csv').query('baseline_scatter_cut==10').copy()
r['east']=(r.registered_ra_deg-150.0634341)*np.cos(np.deg2rad(2.3552216))*3.6e6;r['north']=(r.registered_dec_deg-2.3552216)*3.6e6
q['east']=q.registered_east_mas;q['north']=q.registered_north_mas;q['program']=7814
n=pd.concat([r,q]).merge(f,on='mid_mjd',validate='one_to_one').groupby('program')[['mid_mjd','px','py','east','north']].median().sort_values('mid_mjd')
s=pd.DataFrame(js('hst_refined_fit.zip','summary.json'));e=pd.DataFrame(js('hst_refined_fit.zip','exposures.json'))
ra,dec=np.deg2rad([150.0634341,2.3552216]);E=np.array([-np.sin(ra),np.cos(ra),0]);N=np.array([-np.cos(ra)*np.sin(dec),-np.sin(ra)*np.sin(dec),np.cos(dec)])
if (root/'hst_factors.csv').exists():
 fac=pd.read_csv(root/'hst_factors.csv');mode='supplied exposure factors'
else:
 import swisseph as swe
 times=e.mid_mjd.unique();xyz=np.array([swe.calc_ut(t+2400000.5,swe.EARTH,swe.FLG_MOSEPH|swe.FLG_HELCTR|swe.FLG_J2000|swe.FLG_EQUATORIAL|swe.FLG_XYZ|swe.FLG_TRUEPOS)[0][:3] for t in times]);fac=pd.DataFrame(dict(mid_mjd=times,px=-xyz@E,py=-xyz@N));mode='approximate heliocentric Earth (not HST barycentric)'
H=[];Cs=[];shift=[]
for date in ['2011-12-10','2012-01-26']:
 a=e.query('visit==@date and subset=="common3" and width==1');b=s.query('visit==@date and subset=="common3" and width==1');c=s.query('visit==@date and subset=="common4" and width==1')
 t=a.mid_mjd.mean();dt=(t-60308.94)/365.25;prior=np.array([2.57-70.53*dt,17.59-61.06*dt]);v=b[['dx_mas','dy_mas']].mean().to_numpy();d=a.merge(fac,on='mid_mjd',validate='many_to_one')[['px','py']].mean().to_numpy();H.append([t,*d,*(prior+v)])
 Cs.append(np.cov(a[['dx_mas','dy_mas']].to_numpy(),rowvar=False)/len(a));shift.extend(v-c[['dx_mas','dy_mas']].mean().to_numpy())
rows=np.vstack([H,n.to_numpy()]);labels=['HST2011','HST2012','PRIMER','COSMOS-Web','COSMOS-3D','F140M'];y=rows[:,3:].ravel();t=(rows[:,0]-60308.94)/365.25;A=np.zeros((12,5));A[0::2,0]=1;A[1::2,1]=1;A[0::2,2]=t;A[1::2,3]=t;A[:,4]=rows[:,1:3].ravel();U=np.zeros((12,2));U[:4]=np.tile(np.eye(2),(2,1));v=np.r_[shift,np.zeros(8)]
def fit(a,y,c):
 L=np.linalg.cholesky(c);aw=np.linalg.solve(L,a);yw=np.linalg.solve(L,y);p=np.linalg.lstsq(aw,yw,rcond=None)[0];V=np.linalg.inv(aw.T@aw);chi=np.sum((yw-aw@p)**2);return p,V,float(chi)
out=[]
for hn,nn,fn,shared in [(5,3,8,0),(10,5,10,10),(20,5,15,20),(30,8,20,40)]:
 C=np.diag(np.repeat([hn,hn,nn,nn,nn,fn],2)**2).astype(float);C[:2,:2]+=Cs[0];C[2:4,2:4]+=Cs[1];C+=np.outer(v,v)+shared**2*(U@U.T)
 p,V,ch=fit(A,y,C);_,_,ch0=fit(A[:,:4],y,C);loo={}
 for i,l in enumerate(labels):
  k=np.delete(np.arange(12),[2*i,2*i+1]);pp,VV,_=fit(A[k],y[k],C[np.ix_(k,k)]);loo[l]=[pp[4],np.sqrt(VV[4,4])]
 out.append(dict(floors=[hn,nn,fn,shared],parallax=float(p[4]),sigma=float(np.sqrt(V[4,4])),pm=p[2:4].tolist(),chi2=ch,delta_chi2=ch0-ch,loo=loo))
C=np.diag(np.repeat([10,10,5,5,5,10],2)**2).astype(float);C[:2,:2]+=Cs[0];C[2:4,2:4]+=Cs[1];C+=np.outer(v,v);p,V,ch=fit(np.c_[A,U],y,C)
result=dict(factor_mode=mode,labels=labels,rows=rows.tolist(),hst_reference_shift=shift,rank=int(np.linalg.matrix_rank(A)),scenarios=out,free_HST_offset=dict(parallax=float(p[4]),sigma=float(np.sqrt(V[4,4])),chi2=ch),validated=False)
print(json.dumps(result,indent=2));(root/'joint_pilot_result.json').write_text(json.dumps(result,indent=2))
