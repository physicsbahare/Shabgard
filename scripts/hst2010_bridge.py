"""Expand HST 2010 registration using 101 actual NIRCam anchor candidates.
Use 2011/12 registered HST as a same-instrument bridge, not independent copies.
Outputs diagnostic positions only; no stationarity claim.
"""
import csv,json
from collections import defaultdict
from pathlib import Path
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
import hst_wfc3_ir_epoch as h
OUT=Path('hst2010_bridge');OUT.mkdir(exist_ok=True)
def read(root,name):
 p=list(Path(root).rglob(name));assert len(p)==1,p
 return list(csv.DictReader(open(p[0])))
def dump(name,x):
 (OUT/name).write_text(json.dumps(x,indent=2,allow_nan=False))
def main():
 anchors=h.read('anchor_refs.csv');regs={q['filename']:q for q in h.read('registration_pilot.csv')};stable=h.refs()
 pos={}
 for q in anchors:
  k=q['id']
  if k in stable:pos[k]=stable[k];continue
  if q['anchor_filename'] not in regs:continue
  r=regs[q['anchor_filename']]
  pos[k]=SkyCoord(float(q['ra']),float(q['dec']),unit='deg').spherical_offsets_by(-float(r['east_mas'])*u.mas,-float(r['north_mas'])*u.mas)
 print('EXPANDED_CANDIDATES',len(pos),flush=True)
 rows=[q for q in read('pilot','exposure_audit.csv') if q['status']=='pixel_likelihood']
 old=defaultdict(dict)
 for q in read('pilot','references.csv'):old[q['filename']][q['ref_id']]=np.array([float(q['east_mas']),float(q['north_mas'])])
 measured=[];files={}
 for q in rows:
  n=q['filename'];p=h.download(n,'mast:HST/product/'+n);files[n]=p
  with fits.open(p,memmap=False) as hd:
   sh=hd['SCI',1];w=WCS(sh.header,hd)
   c=np.median([old[n][k] for k in (['6','10'] if q['visit']=='2010-11-03' else ['6','10','11'])],axis=0)
   for k,coord in pos.items():
    try:
     shifted=coord.spherical_offsets_by(c[0]*u.mas,c[1]*u.mas)
     d,s=h.ref_fit(sh.data,hd['ERR',1].data,hd['DQ',1].data,w,shifted)
     v=h.offset(coord,shifted.spherical_offsets_by(d[0]*u.mas,d[1]*u.mas))
     measured.append(dict(filename=n,visit=q['visit'],filter=q['filter'],ref=k,stable_nircam=k in stable,snr=s,east=float(v[0]),north=float(v[1]),residual_east=float(v[0]-c[0]),residual_north=float(v[1]-c[1])))
    except Exception:continue
  dump('measurements.json',measured)
  print('BRIDGE_EXPOSURE',n,sum(x['filename']==n for x in measured),flush=True)
 # Reference positions inferred using only 2011/12, with each visit weighted equally.
 by=defaultdict(lambda:defaultdict(list))
 for q in measured:
  if q['visit']!='2010-11-03':by[q['ref']][q['visit']].append([q['residual_east'],q['residual_north']])
 bridges={}
 for k,d in by.items():
  if len(d)!=2 or min(map(len,d.values()))<2:continue
  z=np.array([np.median(v,axis=0) for v in d.values()]);delta=float(np.linalg.norm(z[0]-z[1]))
  if delta>40:continue
  bridges[k]=dict(residual=np.mean(z,axis=0).tolist(),visit_difference=delta)
 dump('bridge_refs.json',bridges)
 gridaxis=np.arange(-.25,.251,.005);xx,yy=np.meshgrid(gridaxis,gridaxis);grid=np.c_[xx.ravel(),yy.ravel()]
 output=[]
 for limit in [20,40]:
  common=set(k for k,b in bridges.items() if b['visit_difference']<=limit)
  visits=[q for q in rows if q['visit']=='2010-11-03']
  common.intersection_update(*[set(v['ref'] for v in measured if v['filename']==q['filename']) for q in visits])
  print('COMMON_2010',limit,sorted(common),flush=True)
  if len(common)<3:continue
  for excluded in [None,*sorted(common)]:
   ids=sorted(common-{excluded})
   if len(ids)<3:continue
   like=[];rec=[]
   for q in visits:
    n=q['filename'];vv={v['ref']:np.array([v['east'],v['north']])-bridges[v['ref']]['residual'] for v in measured if v['filename']==n and v['ref'] in ids}
    z=np.array([vv[k] for k in ids]);c=np.median(z,axis=0)
    if np.max(np.linalg.norm(z-c,axis=1))>100:continue
    with fits.open(files[n],memmap=False) as hd:
     sh=hd['SCI',1];w=WCS(sh.header,hd)
     center=h.pred(float(q['mid_mjd'])).spherical_offsets_by(c[0]*u.mas,c[1]*u.mas)
     f,v=h.scan(sh.data,hd['ERR',1].data,hd['DQ',1].data,w,center,grid,.60)
     sn=h.combined([(f,v)]);i=np.nanargmax(sn);like.append((f,v))
     rec.append(dict(filename=n,mid_mjd=float(q['mid_mjd']),translation=c.tolist(),ref_scatter=(1.4826*np.median(abs(z-c),axis=0)).tolist(),dx_mas=(1000*grid[i]).tolist(),snr=float(sn[i])))
   if len(like)!=len(visits):continue
   sn=h.combined(like);i=np.nanargmax(sn);o=dict(limit=limit,excluded=excluded,refs=ids,exposures=rec,peak_offset_mas=(1000*grid[i]).tolist(),snr=float(sn[i]),boundary=bool(np.max(abs(grid[i]))>.245),validated=False)
   output.append(o);print('BRIDGE_RESULT',json.dumps(o),flush=True)
   if excluded is None:np.savez_compressed(OUT/f'likelihood_{limit}.npz',grid=grid,snr=sn)
   dump('registration.json',output)
 dump('registration.json',output)
 if not output:print('No viable common 2010 reference set; inspect measured rows',flush=True)
if __name__=='__main__':main()
