"""Run from workspace root: python analysis_outputs/manuscript_sample_date_ordination/make_ordination.py"""
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import json, re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.spatial.distance import pdist, squareform
from threadpoolctl import threadpool_limits
threadpool_limits(limits=1)
warnings.filterwarnings('ignore', module='rdata')
OUT=Path('analysis_outputs/nt_sample_date')
OUT.mkdir(parents=True,exist_ok=True)
BOOK=Path('CFARM_Co_2022_sample_control_list_2023.05.12.xlsx')
b=pd.read_excel(BOOK,sheet_name='sample control list')
b=b[b.UID_code.astype(str).str.fullmatch(r'CS\.\d+')].copy()
assert not b.UID_code.duplicated().any()
b['sample_date']=pd.to_datetime(b.date_sampled,format='%Y.%m.%d',errors='raise').dt.strftime('%Y-%m-%d')
allmeta={}; counts={}; audit=[]
for domain,file in [('16S','combined_16S_ps.rds'),('ITS','combined_its.rds'),('AMF','combined_amf.rds')]:
 info=json.loads(Path(f'manuscript/audit/{domain}_matrix_metadata.json').read_text())
 v=np.load(f'manuscript/audit/{domain}_counts.npz')['values'].astype(float)
 assert not info['taxa_are_rows'][0]
 matrix_ids=info['dimnames'][0]
 m=pd.read_csv(f'manuscript/audit/{domain}_metadata.csv')
 assert v.shape==tuple(info['dim']) and len(m)==len(matrix_ids)
 m['sample_id']=matrix_ids
 uid='UID_code' if 'UID_code' in m else 'Sample_UID_code'
 assert all(str(u).lower().replace('.', '')==re.match(r'cs[0-9]+', sid.lower()).group() for u,sid in zip(m[uid],matrix_ids))
 assert all('_'+str(f).lower()+'_' in sid.lower() for f,sid in zip(m.treatment,matrix_ids))
 m=m.rename(columns={uid:'soil_uid'}).merge(b[['UID_code','sample_date','Site_Name','land_use','Depth_cm']],left_on='soil_uid',right_on='UID_code',how='left',validate='many_to_one',sort=False)
 assert list(m.sample_id)==matrix_ids
 assert np.isfinite(v).all() and (v>=0).all()
 selected=m.treatment.str.lower().eq('mf') & m.management.str.lower().isin(['nt','vtill'])
 expected=pd.read_csv(f'manuscript/audit/{domain}_sample_alpha.csv').iloc[:,0]
 assert set(m.loc[selected,'sample_id'])==set(expected), f'{domain}: manuscript sample mismatch'
 summary=pd.read_csv('manuscript/supplementary_data/dataset_summary.csv').set_index('domain').loc[domain]
 assert selected.sum()==summary.n_samples
 assert (selected & m.management.str.lower().eq('nt')).sum()==summary.n_nt
 assert (selected & m.management.str.lower().eq('vtill')).sum()==summary.n_till
 m['management_group']=m.management.str.lower().map({'nt':'NT','vtill':'till'})
 keep=selected.to_numpy() & m.sample_date.notna().to_numpy() & (v.sum(1)>0)
 assert keep.sum()==selected.sum(), 'A manuscript sample is missing date/counts' 
 audit.append(dict(marker=domain,input_samples=len(m),missing_date=int(m.sample_date.isna().sum()),zero_count_samples=int((v.sum(1)==0).sum()),included=int(keep.sum())))
 m.loc[~keep].to_csv(OUT/f'{domain}_outside_manuscript_cohort.csv',index=False)
 m=m.loc[keep].reset_index(drop=True);v=v[keep];v=v/v.sum(1)[:,None]
 allmeta[domain]=m;counts[domain]=v

from scipy.stats import f_oneway
from statsmodels.stats.multitest import multipletests
B=9999; SEED=20260929
rows=[]; manifests=[]
def permanova(d2,g):
 n=len(g);levels=np.unique(g);total=d2.sum()/(2*n)
 within=sum(d2[np.ix_(g==v,g==v)].sum()/(2*np.sum(g==v)) for v in levels)
 between=total-within
 return (between/(len(levels)-1))/(within/(n-len(levels))),between/total

def dispersion(X,g):
 groups=[X[g==v] for v in np.unique(g)]
 z=[np.linalg.norm(x-x.mean(0),axis=1) for x in groups]
 means=np.array([a.mean() for a in z]);ns=np.array([len(a) for a in z]);grand=np.sum(means*ns)/ns.sum()
 between=np.sum(ns*(means-grand)**2);within=sum(np.sum((a-a.mean())**2) for a in z)
 return (between/(len(z)-1))/(within/(ns.sum()-len(z)))


for scope in ['NT_all','NT_Co_NT_only']:
 for k,domain in enumerate(allmeta):
  mask=allmeta[domain].management_group.eq('NT')
  if scope=='NT_Co_NT_only': mask &= allmeta[domain].site.eq('Co_NT')
  m=allmeta[domain].loc[mask].reset_index(drop=True);v=counts[domain][mask];
  manifests.append(m.assign(marker=domain,scope=scope))
  d=squareform(pdist(v,metric='braycurtis'));d2=d*d;g=pd.Categorical(m.sample_date,categories=sorted(m.sample_date.unique())).codes
  A=-.5*d2;G=A-A.mean(0)-A.mean(1)[:,None]+A.mean();ev=np.linalg.eigvalsh(G);c=max(0,-ev[0]);dc2=d2+2*c;np.fill_diagonal(dc2,0);A=-.5*dc2;G=A-A.mean(0)-A.mean(1)[:,None]+A.mean();ev,U=np.linalg.eigh(G);X=U[:,ev>1e-10]*np.sqrt(ev[ev>1e-10])
  assert np.allclose(squareform(pdist(X))**2,dc2,atol=1e-9)
  F,R2=permanova(d2,g);DF=dispersion(X,g)
  # Independent verification of pseudo-F by projection matrices.
  Z=np.eye(len(np.unique(g)))[g];H=Z@np.linalg.inv(Z.T@Z)@Z.T;J=np.ones_like(H)/len(g)
  A=-.5*d2;rawG=A-A.mean(0)-A.mean(1)[:,None]+A.mean()
  alt=np.trace((H-J)@rawG)/(Z.shape[1]-1)/(np.trace((np.eye(len(g))-H)@rawG)/(len(g)-Z.shape[1]))
  assert np.isclose(F,alt)
  assert np.isclose(DF,f_oneway(*[np.linalg.norm(X[g==v]-X[g==v].mean(0),axis=1) for v in np.unique(g)]).statistic)
  profiles=sorted(m.Site_Name.unique());pindex={p:i for i,p in enumerate(profiles)};ri=np.array([pindex[p] for p in m.Site_Name]);pg=[];sig=[]
  for p in profiles:
   z=m[m.Site_Name==p];assert z.sample_date.nunique()==1
   assert not z.Depth_cm.duplicated().any()
   pg.append(g[z.index[0]]);sig.append('|'.join(sorted(z.Depth_cm)))
  pg=np.array(pg);sig=np.array(sig);strata=[np.where(sig==s)[0] for s in np.unique(sig)]
  rng=np.random.default_rng(SEED+k);exceed=np.zeros(2,int)
  for _ in range(B):
   shuffled=pg.copy()
   for ids in strata:shuffled[ids]=rng.permutation(pg[ids])
   perm=shuffled[ri]
   assert np.array_equal(np.bincount(perm),np.bincount(g))
   exceed += [permanova(d2,perm)[0]>=F-1e-12,dispersion(X,perm)>=DF-1e-12]
  pv=(exceed+1)/(B+1)
  rows.append(dict(scope=scope,marker=domain,n=len(m),profiles=len(profiles),date_groups=len(np.unique(g)),df_date=len(np.unique(g))-1,df_residual=len(m)-len(np.unique(g)),pseudo_F=F,R2=R2,permanova_p=pv[0],dispersion_F=DF,dispersion_p=pv[1],permutations=B,seed=SEED+k))
  print(domain,rows[-1],flush=True)
stats=pd.DataFrame(rows)
for scope,ix in stats.groupby('scope').groups.items():
 for test in ['permanova','dispersion']:
  stats.loc[ix,test+'_q_BH']=multipletests(stats.loc[ix,test+'_p'],method='fdr_bh')[1]
stats.to_csv(OUT/'date_group_statistics.csv',index=False)
manifest=pd.concat(manifests,ignore_index=True)
manifest.to_csv(OUT/'sample_manifest.csv',index=False)
manifest.groupby(['scope','marker','sample_date','site']).agg(n=('sample_id','size'),profiles=('Site_Name','nunique')).reset_index().to_csv(OUT/'sample_counts.csv',index=False)
print(stats.to_string(index=False))
print(manifest.groupby(['scope','marker','sample_date','site']).size().to_string())
