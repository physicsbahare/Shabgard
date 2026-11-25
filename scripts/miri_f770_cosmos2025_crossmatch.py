"""Retrieve public COSMOS2025 local photometry and enumerate all plausible NIRCam counterparts.
No automatic declaration of stationary galaxies or a cross-epoch astrometric tie.
"""
import argparse,csv,json
from pathlib import Path
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.vizier import Vizier

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('refs');ap.add_argument('out');args=ap.parse_args()
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    refs=list(csv.DictReader(open(args.refs)))
    center=SkyCoord(150.0634341,2.3552216,unit='deg')
    radius=max(center.separation(SkyCoord(float(r['ra_deg']),float(r['dec_deg']),unit='deg')).arcsec for r in refs)+3.
    failures=[];table=None
    for server in ['vizier.cds.unistra.fr','vizier.cfa.harvard.edu']:
        try:
            v=Vizier(columns=['**','_RAJ2000','_DEJ2000'],row_limit=-1,vizier_server=server,timeout=120)
            tables=v.query_region(center,radius=radius*u.arcsec,catalog='J/A+A/704/A339/phot',cache=False)
            if len(tables)!=1 or len(tables[0])==0: raise ValueError(f'expected one nonempty phot table, got {len(tables)}')
            table=tables[0];break
        except Exception as e: failures.append({'server':server,'error':repr(e)});print(failures[-1],flush=True)
    if table is None: raise RuntimeError(f'COSMOS2025 public query failed: {failures}')
    table.write(out/'cosmos2025_local.ecsv',format='ascii.ecsv',overwrite=True)
    pairs=[('RAJ2000','DEJ2000'),('_RAJ2000','_DEJ2000'),('RAdeg','DEdeg'),('ra','dec'),('RA','DEC')]
    ra,dec=next(((r,d) for r,d in pairs if r in table.colnames and d in table.colnames),(None,None))
    if ra is None: raise ValueError(f'unknown coordinate columns: {table.colnames}')
    def angle(k):
        col=table[k]
        if col.unit is None:
            # Numeric degrees only; refuse sexagesimal guesses.
            return np.asarray(col,dtype=float)*u.deg
        return u.Quantity(col).to(u.deg)
    cat=SkyCoord(angle(ra),angle(dec))
    finite=np.isfinite(cat.ra.deg)&np.isfinite(cat.dec.deg)
    idcol=next((k for k in ['Id','ID','id','ID_CW','ID_COSMOS','IDcosmos','IDcosmos2025'] if k in table.colnames),None)
    matches=[];counts=[]
    for ref in refs:
        p=SkyCoord(float(ref['ra_deg']),float(ref['dec_deg']),unit='deg')
        sep=p.separation(cat).arcsec
        order=np.argsort(np.where(finite,sep,np.inf))
        cand=[int(i) for i in order if finite[i] and sep[i]<=1.0]
        counts.append(dict(source_id=ref['source_id'],n_within_0p3=sum(sep[i]<=.3 for i in cand),
            n_within_0p5=sum(sep[i]<=.5 for i in cand),n_within_1p0=len(cand),
            nearest_arcsec=float(sep[order[0]]),second_arcsec=float(sep[order[1]]) if len(order)>1 else None))
        for rank,i in enumerate(cand,1):
            dra,ddec=p.spherical_offsets_to(cat[i])
            matches.append(dict(source_id=ref['source_id'],catalog_row=i,catalog_id=str(table[idcol][i]) if idcol else '',
                rank=rank,separation_arcsec=float(sep[i]),delta_ra_cosdec_mas=dra.mas,delta_dec_mas=ddec.mas,
                catalog_ra_deg=cat[i].ra.deg,catalog_dec_deg=cat[i].dec.deg,
                classification='unvalidated_counterpart'))
    def write(name,rows):
        if rows:
            with (out/name).open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    write('counterpart_candidates.csv',matches);write('counterpart_counts.csv',counts)
    summary=dict(query_catalog='J/A+A/704/A339/phot',server=server,query_radius_arcsec=radius,n_catalog=len(table),columns=table.colnames,
        id_column=idcol,coordinate_columns=[ra,dec],n_references=len(refs),n_within_0p3=sum(r['n_within_0p3']>0 for r in counts),
        n_unique_within_0p3=sum(r['n_within_0p3']==1 for r in counts),n_within_0p5=sum(r['n_within_0p5']>0 for r in counts),
        n_ambiguous_within_0p5=sum(r['n_within_0p5']>1 for r in counts),n_within_1p0=sum(r['n_within_1p0']>0 for r in counts),
        failed_attempts=failures,stationary_galaxy_classification_validated=False,cross_epoch_registration_validated=False)
    (out/'summary.json').write_text(json.dumps(summary,indent=2,default=lambda x: x.item() if isinstance(x,np.generic) else str(x)));print(json.dumps(summary,indent=2,default=lambda x: x.item() if isinstance(x,np.generic) else str(x)),flush=True)
if __name__=='__main__':main()
