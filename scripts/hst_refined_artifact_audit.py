import zipfile,json,sys
import numpy as np
z=zipfile.ZipFile(sys.argv[1])
s=json.load(z.open('summary.json')); e=json.load(z.open('exposures.json'))
for date in ['2011-12-10','2012-01-26']:
 def get(f,k,w=1):return next(r for r in s if r['visit']==date and r['filter']==f and r['subset']==k and r['width']==w)
 def v(r):return np.array([r['dx_mas'],r['dy_mas']])
 for k in ['common4','common3']:
  a,b=[get(f,k) for f in ['F125W','F160W']]; delta=v(a)-v(b); C=np.array(a['formal_cov_mas2'])+np.array(b['formal_cov_mas2']);print(date,k,'filter_delta',delta,'norm',np.linalg.norm(delta),'formal_only_sigma',np.sqrt(delta@np.linalg.solve(C,delta)))
  for f in ['F125W','F160W']:
   q=np.array([v(r) for r in e if r['visit']==date and r['filter']==f and r['subset']==k and r['width']==1]);print(f,'dither_difference',q[1]-q[0],'width_span',np.ptp([v(get(f,k,w)) for w in [.85,1,1.15]],axis=0))
 print('subset_shift_common3_minus4',np.mean([v(get(f,'common3'))-v(get(f,'common4')) for f in ['F125W','F160W']],axis=0))
