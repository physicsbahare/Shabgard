Artifact 15169064227 downloaded and analyzed with Python. Seven accepted references, four exposures each. For each reference subtract its own temporal mean (removes static offsets, not an absolute-frame test). Predict exposure common mode using median of the other six centered references, subtract from held-out reference. Pooled RMS east/north before 12.44981133/8.78634159 mas; held-out after 10.94106310/7.85688566 mas (12.1%/10.6% reduction). Individual RMS before -> after, east/north: 10 (4.86,9.89)->(11.32,17.85); 36 (8.36,8.24)->(3.65,2.03); 41 (13.62,8.11)->(10.07,1.60); 43 (25.11,6.05)->(21.75,5.55); 6 (4.92,9.41)->(4.82,2.47); 64 (10.76,8.23)->(8.12,4.17); 92 (5.91,10.78)->(5.76,7.25). Ref 10 worsens markedly. Selection and temporal centering use same visit, so not fully independent external validation. No target correction adopted. Raw FLT pixels absent in artifact; physical injection remains pending. Pilot unchanged.

Reproduction (extract artifact first):
```python
import json,pandas as pd,numpy as np
d=pd.DataFrame(json.load(open('measurements.json')))
r=json.load(open('registration.json'))[0]
a=d[(d.visit=='2010-11-03')&d.ref.isin(r['refs'])]
x=np.stack([a.pivot(index='filename',columns='ref',values=k).values for k in ['residual_east','residual_north']],axis=-1)
x-=x.mean(axis=0)
y=np.stack([x[:,i]-np.median(np.delete(x,i,axis=1),axis=1) for i in range(x.shape[1])],axis=1)
print(np.sqrt(np.mean(x*x,axis=(0,1))),np.sqrt(np.mean(y*y,axis=(0,1))))
print(np.sqrt(np.mean(x*x,axis=0)),np.sqrt(np.mean(y*y,axis=0)))
```
