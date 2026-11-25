"""Relative same-visit translation network on real F770W field candidates.
No absolute reference, morphology classification, or cross-epoch tie.
"""
import csv,json,itertools,sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

def analyze(path,radius=500.,snr=0.,nboot=60,export=None):
    rows=[r for r in csv.DictReader(open(path)) if float(r['fit_snr_proxy'])>=snr]
    files=sorted(set(r['filename'] for r in rows)); idx=np.array([files.index(r['filename']) for r in rows]); nf=len(files)
    xy=np.array([[(float(r['ra'])-150.0634341)*np.cos(np.deg2rad(2.3552216))*3.6e6,(float(r['dec'])-2.3552216)*3.6e6] for r in rows])
    parent=list(range(len(rows)))
    def root(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]]; i=parent[i]
        return i
    edges=[]
    for i,j in itertools.combinations(range(nf),2):
        a,b=np.where(idx==i)[0],np.where(idx==j)[0]
        d,k=cKDTree(xy[b]).query(xy[a]); _,l=cKDTree(xy[a]).query(xy[b])
        p=[(a[u],b[v]) for u,v in enumerate(k) if d[u]<radius and l[v]==u]
        if p: edges.append((i,j,len(p)))
        for u,v in p: parent[root(u)]=root(v)
    groups={}
    for i in range(len(rows)): groups.setdefault(root(i),[]).append(i)
    bad=[g for g in groups.values() if len(set(idx[g]))!=len(g)]
    gs=[np.array(g) for g in groups.values() if len(g)>1 and len(set(idx[g]))==len(g)]
    ref=max(range(nf),key=lambda i:sum(i in idx[g] for g in gs))
    E=np.eye(nf)[idx][:,np.arange(nf)!=ref]
    def solve(groups):
        weights=[np.ones(len(g)) for g in groups]
        for _ in range(5):
            A=[];B=[]
            for g,w in zip(groups,weights):
                z=E[g];y=xy[g];A.append((z-np.average(z,axis=0,weights=w))*np.sqrt(w[:,None]));B.append((y-np.average(y,axis=0,weights=w))*np.sqrt(w[:,None]))
            A=np.concatenate(A);B=np.concatenate(B)
            if np.linalg.matrix_rank(A)<nf-1: raise ValueError('disconnected graph')
            b=np.linalg.lstsq(A,B,rcond=None)[0]
            residual=[]
            for g,w in zip(groups,weights):
                z=xy[g]-E[g]@b;residual.append(z-np.average(z,axis=0,weights=w))
            norms=np.concatenate([np.linalg.norm(r,axis=1) for r in residual]); scale=max(10.,float(np.median(norms))/1.1774)
            weights=[np.minimum(1.,1.5*scale/np.maximum(np.linalg.norm(r,axis=1),1e-12)) for r in residual]
        offset=np.insert(b,ref,0.,axis=0)
        return offset,np.concatenate(residual)
    off,res=solve(gs);rng=np.random.default_rng(639188);boot=[]
    for _ in range(nboot):
        try: boot.append(solve([gs[k] for k in rng.integers(len(gs),size=len(gs))])[0])
        except ValueError: pass
    loo=[]
    for k in range(len(gs)):
        try: loo.append(solve(gs[:k]+gs[k+1:])[0])
        except ValueError: pass
    if export is not None:
        export=Path(export);export.mkdir(parents=True,exist_ok=True)
        def write(name,records):
            with (export/name).open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
        refs=[];members=[]
        for j,g in enumerate(gs):
            z=xy[g]-off[idx[g]]
            center=np.median(z,axis=0)
            sid=f'F770REF{j+1:03d}'
            ra=150.0634341+center[0]/(3.6e6*np.cos(np.deg2rad(2.3552216)))
            dec=2.3552216+center[1]/3.6e6
            refs.append(dict(source_id=sid,ra_deg=ra,dec_deg=dec,n_dithers=len(g),
                median_snr_proxy=float(np.median([float(rows[k]['fit_snr_proxy']) for k in g])),
                median_profile_ratio=float(np.median([float(rows[k]['profile_ratio']) for k in g])),
                classification='unclassified',frame='relative_F770W_visit'))
            for k,p in zip(g,z):
                members.append(dict(source_id=sid,filename=rows[k]['filename'],x=rows[k]['x'],y=rows[k]['y'],
                    ra_raw_deg=rows[k]['ra'],dec_raw_deg=rows[k]['dec'],ra_residual_mas=p[0]-center[0],dec_residual_mas=p[1]-center[1]))
        write('reference_sources.csv',refs);write('reference_members.csv',members)
        if len(boot)>1:
            samples=np.array(boot).reshape(len(boot),-1)
            np.savez(export/'relative_translation_covariance.npz',covariance_mas2=np.cov(samples,rowvar=False),
                samples_mas=samples,filenames=np.array(files),axis_order=np.array(['RA_cosdec','Dec']),reference=files[ref])
    mad=1.4826*np.median(np.abs(res-np.median(res,axis=0)),axis=0)
    return dict(radius_mas=radius,snr_cut=snr,reference=files[ref],n_files=nf,n_edges=len(edges),pair_links=sum(e[2] for e in edges),n_conflicting_groups=len(bad),n_sources=len(gs),n_observations=sum(map(len,gs)),residual_mad_mas=mad.tolist(),bootstrap_valid=len(boot),bootstrap_attempts=nboot,loo_valid=len(loo),edges=[dict(a=files[i],b=files[j],n=n) for i,j,n in edges],offsets=[dict(filename=f,n_sources=sum(i in idx[g] for g in gs),ra_mas=off[i,0],dec_mas=off[i,1],bootstrap_sd_mas=np.std(boot,axis=0,ddof=1)[i].tolist() if len(boot)>1 else None,loo_max_shift_mas=float(np.max(np.linalg.norm(np.array(loo)[:,i]-off[i],axis=1))) if loo else None) for i,f in enumerate(files)])

if __name__=='__main__':
    path=sys.argv[1];out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
    results=[]
    for radius,snr in [(500.,0.),(250.,0.),(150.,0.),(250.,15.),(250.,30.)]:
        try: results.append(analyze(path,radius,snr,nboot=300 if radius==250 and snr==0 else 60,export=out if radius==250 and snr==0 else None))
        except ValueError as e: results.append(dict(radius_mas=radius,snr_cut=snr,error=str(e)))
    (out/'network.json').write_text(json.dumps(results,indent=2));print(json.dumps(results,indent=2))
