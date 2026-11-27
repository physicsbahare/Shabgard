"""Forecast additional independent epoch leverage with existing exact 6-epoch fit covariance.
No candidate detection or hypothetical centroid inserted into the data.
Run with joint_exact.zip as first argument; requires numpy pandas pyswisseph.
"""
import json,zipfile,numpy as np,pandas as pd,swisseph as swe,sys
from datetime import datetime,timezone
z=zipfile.ZipFile(sys.argv[1]);j=json.load(z.open('joint_pilot_result.json'));f=pd.read_csv(z.open('hst_factors.csv'));rows=np.array(j['rows'])
ra,dec=np.deg2rad([150.0634341,2.3552216]);E=np.array([-np.sin(ra),np.cos(ra),0]);N=np.array([-np.cos(ra)*np.sin(dec),-np.sin(ra)*np.sin(dec),np.cos(dec)])
def earth(date):
 t=datetime.fromisoformat(date).replace(tzinfo=timezone.utc).timestamp()/86400+40587+.5
 v=np.array(swe.calc_ut(t+2400000.5,swe.EARTH,swe.FLG_MOSEPH|swe.FLG_HELCTR|swe.FLG_J2000|swe.FLG_EQUATORIAL|swe.FLG_XYZ|swe.FLG_TRUEPOS)[0][:3]);return [t,-v@E,-v@N]
def design(x):
 t=(x[:,0]-60308.94)/365.25;A=np.zeros((2*len(x),5));A[::2,0]=1;A[1::2,1]=1;A[::2,2]=t;A[1::2,3]=t;A[:,4]=x[:,1:3].ravel();return A
C=np.diag(np.repeat([10,10,5,5,5,10],2)**2).astype(float);C[:2,:2]+=[[43.2292,-7.2917],[-7.2917,10.4167]];C[2:4,2:4]+=[[30.7292,25],[25,20.8333]];v=np.r_[j['hst_reference_shift'],np.zeros(8)];C+=np.outer(v,v);U=np.zeros((12,2));U[:4]=np.tile(np.eye(2),(2,1));C+=100*U@U.T
A=design(rows);base=np.sqrt(np.linalg.inv(A.T@np.linalg.solve(C,A))[4,4]);print('base',base)
opts={'HST2010':(f.loc[f.mid_mjd<55600,['mid_mjd','px','py']].mean().tolist(),True),'HST2021':(earth('2021-11-27'),True),'MIRI2025':(earth('2025-12-29'),False),'NIRCam future winter illustration':(earth('2026-12-29'),False)}
for k,(x,shared) in opts.items():
 B=design(np.vstack([rows[:,:3],x]));print(k,x)
 for s in [2,3,5,10,20,40]:
  D=np.zeros((14,14));D[:12,:12]=C;D[12:,12:]=s*s*np.eye(2)
  if shared:D[12:,12:]+=100*np.eye(2);D[:12,12:]=100*U;D[12:,:12]=100*U.T
  V=np.linalg.inv(B.T@np.linalg.solve(D,B));print(s,round(float(np.sqrt(V[4,4])),4))
