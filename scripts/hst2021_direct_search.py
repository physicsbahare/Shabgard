"""Search 2021 HST visit for associated direct images, then test actual FLT pixels.
No grism trace/zero order counted as direct astrometry.
"""
import csv,json
from pathlib import Path
import numpy as np
from astropy.time import Time
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astroquery.mast import Observations
import hst_wfc3_ir_epoch as h
OUT=Path('hst2021_search');OUT.mkdir(exist_ok=True)
def save(n,r):
 (OUT/n).write_text(json.dumps(r,indent=2,allow_nan=False))
def main():
 t=Time('2021-11-27').mjd;p=h.pred(t)
 queries=[('visit',lambda:Observations.query_criteria(obs_collection='HST',obs_id='iehn8r*')),('region',lambda:Observations.query_region(p,radius=2*u.arcmin))]
 obsrows=[];prods={};errors=[]
 for label,query in queries:
  try:
   obs=query()
   if label=='region':obs=obs[(np.asarray(obs['obs_collection'])=='HST')&(np.asarray(obs['proposal_id']).astype(str)=='16443')]
   for q in obs:
    obsrows.append({k:str(q[k]) for k in ['obs_id','obsid','instrument_name','filters','t_min','t_max','s_region','proposal_id','dataRights']})
   for q in Observations.get_product_list(obs):
    n=str(q['productFilename']).lower()
    if n.endswith(('_flt.fits','_asn.fits')):prods[n]=str(q['dataURI'])
  except Exception as e:errors.append(dict(query=label,error=repr(e)))
 save('observations.json',obsrows);save('products.json',prods)
 print('2021_OBSERVATIONS',json.dumps(obsrows),flush=True)
 # ASN may supply direct exposure names; restrict to the relevant visit.
 direct={k:v for k,v in prods.items() if k.startswith('iehn8r')}
 products,manifest=h.discover('iehn8r',['iehn8r010'])
 direct.update(products)
 for n,uri in list(direct.items()):
  if not n.endswith('_asn.fits'):continue
  try:
   with fits.open(h.download(n,uri),memmap=False) as hd:
    for q in hd[1].data:
     if 'EXP' in str(q['MEMTYPE']).upper():
      m=str(q['MEMNAME']).strip().lower()+'_flt.fits';direct.setdefault(m,'mast:HST/product/'+m)
  except Exception as e:errors.append(dict(file=n,error=repr(e)))
 # Include other program-16443 local direct visit candidates only if metadata suggests relevant dates.
 for o in obsrows:
  try:
   if abs(float(o['t_min'])-t)>3 or 'G141' in o['filters'] or 'WFC3/IR' not in o['instrument_name']:continue
   for n,uri in prods.items():
    if n.startswith(o['obs_id'][:6].lower()):direct[n]=uri
  except ValueError:continue
 save('attempted_products.json',direct);save('errors.json',errors)
 results=[];refs=h.refs();a=np.arange(-.6,.601,.02);xx,yy=np.meshgrid(a,a);grid=np.c_[xx.ravel(),yy.ravel()]
 for n,uri in sorted(direct.items()):
  if not n.endswith('_flt.fits'):continue
  row=dict(filename=n,status='error')
  try:
   with fits.open(h.download(n,uri),memmap=False) as hd:
    ph=hd[0].header;f=str(ph.get('FILTER','')).upper();row['filter']=f
    if ph.get('DETECTOR')!='IR' or not f.startswith('F'):
     row['status']='not_direct_IR';continue
    sh=hd['SCI',1];w=WCS(sh.header,hd);t0=h.exposure_mjd(ph,sh.header);pos=h.pred(t0)
    x,y=h.pixel(w,pos);row.update(mid_mjd=t0,x=float(x),y=float(y))
    if not (12<x<sh.data.shape[1]-13 and 12<y<sh.data.shape[0]-13):
     row['status']='outside_detector';continue
    ref=[]
    for k,p0 in refs.items():
     try:
      v,s=h.ref_fit(sh.data,hd['ERR',1].data,hd['DQ',1].data,w,p0);ref.append([k,*map(float,v),s])
     except Exception:continue
    row['refs']=ref
    c=np.median(np.array([v[1:3] for v in ref]),axis=0) if len(ref)>=4 else np.zeros(2)
    row['registered_pilot']=len(ref)>=4;row['translation']=c.tolist()
    pos=pos.spherical_offsets_by(c[0]*u.mas,c[1]*u.mas)
    sig=.65 if f in ['F160W','F167N'] else .55
    flux,var=h.scan(sh.data,hd['ERR',1].data,hd['DQ',1].data,w,pos,grid,sig)
    sn=h.combined([(flux,var)]);i=np.nanargmax(sn)
    row.update(status='pixel_likelihood',peak_snr=float(sn[i]),offset_mas=(1000*grid[i]).tolist(),boundary=bool(np.max(abs(grid[i]))>.58))
    np.savez_compressed(OUT/(n+'.npz'),grid=grid,snr=sn)
  except Exception as e:row['error']=repr(e)
  finally:
   results.append(row);save('results.json',results);print('2021_FLT',json.dumps(row),flush=True)
 print('2021_SUMMARY',len(obsrows),len(direct),len(results),flush=True)
if __name__=='__main__':main()
