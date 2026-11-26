"""Screen all 2022 PRIMER F770W CAL headers, then audit actual pixels of plausible members."""
import csv,json
from pathlib import Path
from nircam_same_filter_astrometry import get_header,footprint,download
from select_miri_cal_pilot import archive_products
from miri_calibration_pilot import inspect
OUT=Path('primer2022');OUT.mkdir(exist_ok=True)
def write(name,rows):
    if rows:
        with (OUT/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=sorted(set().union(*(q.keys() for q in rows))));w.writeheader();w.writerows(rows)
def main():
    # Query all public CALs; prior pilot sampled only 3 of 86.
    items={q['filename']:q for q in archive_products('109374093') if q['rights'].upper() in ('PUBLIC','')}
    if not items:raise RuntimeError('No public PRIMER F770W products')
    audit=[];pixels=[];unknown=0
    for name,q in sorted(items.items()):
        uri=q['data_uri'] or 'mast:JWST/product/'+name
        rec=dict(filename=name,status='error')
        try:
            p,s=get_header(uri)
            rec['filter']=p.get('FILTER','')
            if p.get('FILTER')!='F770W':rec['status']='other_filter';continue
            hit=footprint(p,s);rec['header_hit']=str(hit)
            if hit is False:rec['status']='outside_header_polygon';continue
            if hit is None:
                unknown+=1
                if unknown>12:rec['status']='unknown_header_download_cap';continue
            fn=download(uri,name)
            z=inspect(fn,OUT);pixels.append(z);rec['status']=z['status']
            print('PRIMER_PIXEL',json.dumps(z),flush=True)
        except Exception as e:rec['error']=repr(e)
        finally:
            audit.append(rec);write('header_audit.csv',audit);write('pixel_audit.csv',pixels)
    from collections import Counter
    result=dict(total_products=len(items),counts=dict(Counter(q['status'] for q in audit)),pixel_audits=len(pixels),source_association_validated=False)
    print('PRIMER2022_SUMMARY',json.dumps(result),flush=True)
    (OUT/'summary.json').write_text(json.dumps(result,indent=2))
    if all(q['status']=='error' for q in audit):raise RuntimeError('All PRIMER header checks errored')
if __name__=='__main__':main()
