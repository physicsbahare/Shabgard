import json
from pathlib import Path
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
from astroquery.gaia import Gaia
import hst_wfc3_ir_epoch as h
out=Path("gaia_reference_results");out.mkdir(exist_ok=True)
refs=h.refs()
rows=[]
for k,p in refs.items():
    query=f"""SELECT source_id,ra,dec,parallax,parallax_error,pmra,pmdec,pmra_error,pmdec_error,ruwe,phot_g_mean_mag FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{p.ra.deg},{p.dec.deg},0.001))"""
    try:
        tab=Gaia.launch_job(query).get_results()
        matches=[]
        for t in tab:
            c=SkyCoord(t['ra']*u.deg,t['dec']*u.deg)
            sep=c.separation(p).arcsec
            propagated=None
            if not np.ma.is_masked(t['pmra']) and not np.ma.is_masked(t['pmdec']):
                c=SkyCoord(ra=t['ra']*u.deg,dec=t['dec']*u.deg,pm_ra_cosdec=t['pmra']*u.mas/u.yr,pm_dec=t['pmdec']*u.mas/u.yr,obstime=Time(2016.,format='jyear'))
                propagated=c.apply_space_motion(new_obstime=Time(2024.,format='jyear')).separation(p).arcsec
            matches.append(dict(source_id=str(t['source_id']),sep_2016=sep,sep_2024=propagated,**{n:None if np.ma.is_masked(t[n]) else float(t[n]) for n in ['parallax','parallax_error','pmra','pmdec','pmra_error','pmdec_error','ruwe','phot_g_mean_mag']}))
        row=dict(ref_id=k,ra=p.ra.deg,dec=p.dec.deg,matches=matches,status='queried_not_classified')
    except Exception as e:
        row=dict(ref_id=k,status='query_failed',error=repr(e))
    rows.append(row)
    (out/'matches.json').write_text(json.dumps(rows,indent=2))
    print('GAIA_REF',json.dumps(row),flush=True)
if not any(r['status']=='queried_not_classified' for r in rows):raise RuntimeError('No successful Gaia queries')
