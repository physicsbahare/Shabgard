"""Use real archived positions for a differential-only MIRI constraint.
Inputs: the three named ZIPs from runs 45795708087, 45808293268, 45850200442.
No external frame tie, no publication-ready fit. Error scenarios are assumptions.
"""
import argparse,json,zipfile
from pathlib import Path
import numpy as np
import pandas as pd

def main(root):
    def read(z,n):
        with zipfile.ZipFile(root/z) as f:return pd.read_csv(f.open(n))
    r=read('baseline_nircam_new.zip','registration_pilot.csv')
    f=read('fourth_epoch_new.zip','jwst_parallax_factors.csv')
    q=read('fourth_epoch_new.zip','registration.csv').query('baseline_scatter_cut==10').copy()
    m=read('crossvisit_displacement.zip','positions.csv')
    with zipfile.ZipFile(root/'crossvisit_displacement.zip') as z:s=json.load(z.open('summary.json'))
    r['east']=(r.registered_ra_deg-150.0634341)*np.cos(np.deg2rad(2.3552216))*3.6e6
    r['north']=(r.registered_dec_deg-2.3552216)*3.6e6
    r=r.merge(f,on='mid_mjd',validate='one_to_one')
    q=q.merge(f,on='mid_mjd',validate='one_to_one')
    q['east']=q.registered_east_mas;q['north']=q.registered_north_mas;q['program']=7814
    cols=['mid_mjd','px','py','east','north']
    epochs=pd.concat([r,q]).groupby('program')[cols].median().sort_values('mid_mjd')
    t=(epochs.mid_mjd.to_numpy()-60308.94)/365.25
    p=epochs[['px','py']].to_numpy();y=epochs[['east','north']].to_numpy().ravel()
    A=np.zeros((8,5));A[0::2,0]=1;A[1::2,1]=1;A[0::2,2]=t;A[1::2,3]=t;A[:,4]=p.ravel()
    g=m.groupby('visit')[['mid_mjd','px','py']].mean()
    d=g.loc[2024]-g.loc[2022];dt=d.mid_mjd/365.25;dp=d[['px','py']].to_numpy()
    B=np.array([[0,0,dt,0,dp[0]],[0,0,0,dt,dp[1]]])
    z=next(x for x in s['results'] if x['radius_mas']==250 and x['snr_cut']==0)
    obs=np.array(z['displacement_2024_minus_2022_mas'])
    C=np.array(s['translation_only_contrast_covariance_mas2'])+np.diag(np.array(z['scatter_2022_mas'])**2/4+np.array(z['scatter_2024_mas'])**2/8)
    b=np.linalg.lstsq(A[:,:4],y,rcond=None)[0]
    print('epochs',epochs.to_string())
    print('dt dp',dt,dp,'C',C,'sd',np.sqrt(C.diagonal()))
    print('PM only equal weighting',b,'predict',B[:,:4]@b,'resid',obs-B[:,:4]@b)
    print('MIRI PM assuming no parallax',obs/dt)
    for nerr,ferr,merr in [(3,8,0),(3,8,10),(5,10,10),(5,15,20)]:
        D=np.diag(np.repeat([nerr]*3+[ferr],2)**2);G=C+np.eye(2)*2*merr**2
        def calc(a,v,c):
            W=np.linalg.inv(c);V=np.linalg.inv(a.T@W@a);return V@a.T@W@v,V
        b,v=calc(A,y,D)
        full=np.block([[D,np.zeros((8,2))],[np.zeros((2,8)),G]])
        bb,vv=calc(np.vstack([A,B]),np.r_[y,obs],full)
        res=obs-B@b;S=G+B@v@B.T
        print('scenario',nerr,ferr,merr,'pred',B@b,'res',res,'mahal',np.sqrt(res@np.linalg.solve(S,res)),'sigma_pi_ratio',np.sqrt(vv[4,4]/v[4,4]),'rank',np.linalg.matrix_rank(np.vstack([A,B])))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);main(p.parse_args().root)
