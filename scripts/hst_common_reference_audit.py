"""Re-register existing HST coarse-grid exposure peaks using identical reference IDs.
Diagnostic only: grid quantization and shared reference errors remain.
"""
import zipfile,sys
import pandas as pd
import numpy as np
with zipfile.ZipFile(sys.argv[1]) as z:
    a=pd.read_csv(z.open('exposure_audit.csv')).query('registered==True')
    r=pd.read_csv(z.open('references.csv'))
for ids in [[6,10,11,16],[6,10,11],[6,10]]:
    out=[]
    for _,q in a.iterrows():
        v=r[(r.filename==q.filename)&r.ref_id.isin(ids)][['east_mas','north_mas']].to_numpy()
        assert len(v)==len(ids)
        c=np.median(v,axis=0)
        old=q[['translation_east_mas','translation_north_mas']].to_numpy(float)
        d=q[['dx_arcsec','dy_arcsec']].to_numpy(float)*1000+old-c
        out.append([q.visit,q['filter'],*d])
    b=pd.DataFrame(out,columns=['visit','filter','x','y'])
    print(ids,b.groupby(['visit','filter'])[['x','y']].mean().to_string())
    print(b.groupby('visit')[['x','y']].std())
