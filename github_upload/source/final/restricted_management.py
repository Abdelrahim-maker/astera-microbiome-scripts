"""Management-only PERMANOVA; no collection-date input or model."""
from common import ROOT, load_counts, load_metadata
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.linalg import svd
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits
threadpool_limits(limits=1)
B=9999
OUT=ROOT/'reproduced/restricted_permanova'
OUT.mkdir(parents=True,exist_ok=True)
def projector(m, columns):
    X=np.column_stack([pd.get_dummies(m[c],drop_first=True,dtype=float).to_numpy() for c in columns])
    X=X-X.mean(axis=0)
    u,s,_=svd(X,full_matrices=False)
    q=u[:,s>1e-9]
    return q@q.T,q.shape[1]
def gram(d):
    a=-0.5*d*d
    return a-a.mean(0)-a.mean(1)[:,None]+a.mean()
def fit(G,H,k):
    n=len(G); total=np.trace(G); ss=np.sum(H*G)
    r2=ss/total
    return r2,1-(1-r2)*(n-1)/(n-k-1),(ss/k)/((total-ss)/(n-k-1))
def permutations(m, within_site, rng):
    groups={}
    for profile,g in m.groupby('Site_Name'):
        assert g.management_group.nunique()==g.site.nunique()==1
        assert not g.Depth_cm.duplicated().any()
        ids=g.sort_values('Depth_cm').index.to_numpy()
        key=(g.site.iloc[0] if within_site else 'all',tuple(m.loc[ids,'Depth_cm']))
        groups.setdefault(key,[]).append(ids)
    maps=[]
    for _ in range(B):
        ix=np.arange(len(m))
        for group in groups.values():
            for target,source in zip(group,rng.permutation(len(group))):
                ix[target]=group[source]
        assert np.array_equal(np.sort(ix),np.arange(len(m)))
        assert np.array_equal(m.Depth_cm.to_numpy(),m.Depth_cm.to_numpy()[ix])
        if within_site: assert np.array_equal(m.site.to_numpy(),m.site.to_numpy()[ix])
        maps.append(ix)
    return maps

rows=[]
meta=load_metadata()[['marker','sample_id','Site_Name','management_group','site','Depth_cm']]
for i,marker in enumerate(['16S','ITS','AMF']):
 m=meta[meta.marker.eq(marker)].reset_index(drop=True)
 info,counts=load_counts(marker)
 lookup={s:i for i,s in enumerate(info['dimnames'][0])}
 v=counts[[lookup[s] for s in m.sample_id]].astype(float)
 v/=v.sum(1)[:,None]
 G=gram(squareform(pdist(v,metric='braycurtis')))
 H,k=projector(m,['management_group'])
 r2,adj,F=fit(G,H,k)
 seed=20260930+2*i
 maps=permutations(m,True,np.random.default_rng(seed))
 exceed=sum(fit(G[np.ix_(ix,ix)],H,k)[2]>=F-1e-12 for ix in maps)
 rows.append(dict(marker=marker,model='management',method='PERMANOVA_uncorrected_Bray',restriction='site_and_depth_coverage',n=len(m),df=k,R2=r2,adjusted_R2=adj,pseudo_F=F,p=(exceed+1)/(B+1),permutations=B,seed=seed,lingoes_constant=0))
 print(marker,rows[-1]['p'],flush=True)
x=pd.DataFrame(rows);x['q_BH']=multipletests(x.p,method='fdr_bh')[1]
x.to_csv(OUT/'model_statistics.csv',index=False)
print('Completed management-only PERMANOVA; outputs generated locally, no published reference tables required.')
