from common import ROOT, load_counts, load_metadata, load_alpha
"""Restricted whole-profile PERMANOVA and separate management/date dbRDA."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.linalg import svd
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
threadpool_limits(limits=1)
OUT=ROOT/'reproduced/restricted_permanova';OUT.mkdir(parents=True,exist_ok=True)
B=9999
SEED=20260930
m_all=load_metadata()
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
        assert g.management_group.nunique()==g.sample_date.nunique()==g.site.nunique()==1
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
rows=[]; conditional=[]; scores=[]; caches={}; design=[]
for di,marker in enumerate(['16S','ITS','AMF']):
    m=m_all[m_all.marker.eq(marker)].reset_index(drop=True)
    info,counts=load_counts(marker)
    assert not info['taxa_are_rows'][0]
    lookup={s:i for i,s in enumerate(info['dimnames'][0])}
    v=counts[[lookup[s] for s in m.sample_id]].astype(float)
    assert np.all(v.sum(1)>0)
    v/=v.sum(1)[:,None]
    d=squareform(pdist(v,metric='braycurtis')); Graw=gram(d)
    c=max(0.,-np.linalg.eigvalsh(Graw)[0]); d2=d*d+2*c; np.fill_diagonal(d2,0)
    G=gram(np.sqrt(d2)); e,u=np.linalg.eigh(G); ok=e>1e-10; Y=u[:,ok]*np.sqrt(e[ok])
    assert np.allclose(squareform(pdist(Y))**2,d2,atol=1e-9)
    maps={name:permutations(m,flag,np.random.default_rng(SEED+di*2+j)) for j,(name,flag) in enumerate([('site_and_depth_coverage',True),('depth_coverage_only',False)])}
    for name,ixs in maps.items():
        changed=sum(not np.array_equal(m.management_group.to_numpy(),m.management_group.to_numpy()[ix]) for ix in ixs)
        design.append(dict(marker=marker,restriction=name,profiles=m.Site_Name.nunique(),permutations_changing_management=changed,permutations=B))
        assert changed>0
    for label,cols in [('management',['management_group']),('date',['sample_date'])]:
        H,k=projector(m,cols); assert np.allclose(H@H,H)
        for method,g in [('PERMANOVA_uncorrected_Bray',Graw),('dbRDA_Lingoes',G)]:
            r2,adj,F=fit(g,H,k)
            for restriction,ixs in maps.items():
                exceed=sum(fit(g[np.ix_(ix,ix)],H,k)[2]>=F-1e-12 for ix in ixs)
                rows.append(dict(marker=marker,model=label,method=method,restriction=restriction,n=len(m),df=k,R2=r2,adjusted_R2=adj,pseudo_F=F,p=(exceed+1)/(B+1),permutations=B,seed=SEED+di*2+(restriction=='depth_coverage_only'),lingoes_constant=c if method=='dbRDA_Lingoes' else 0))
        fitted=H@Y; _,s,vt=svd(fitted,full_matrices=False)
        assert np.isclose(np.sum(s*s),np.sum(H*G))
        x=Y@vt[0]
        if k==1:
            ur,sr,_=svd(Y-fitted,full_matrices=False); y=ur[:,0]*sr[0]; ylab='Residual PCoA 1'; vy=sr[0]**2
        else:
            y=Y@vt[1]; ylab='dbRDA 2'; vy=s[1]**2
        pct=100*np.array([s[0]**2,vy])/np.trace(G)
        caches[(marker,label)]=(m,x,y,pct,ylab)
        for i,r in m.iterrows():scores.append(dict(marker=marker,model=label,sample_id=r.sample_id,management=r.management_group,date=r.sample_date,x=x[i],y=y[i],y_axis=ylab))
    for focal,other in [('management_group',['sample_date']),('sample_date',['management_group']),('management_group',['sample_date','site','Depth_cm']),('sample_date',['management_group','site','Depth_cm'])]:
        H0,k0=projector(m,other); H1,k1=projector(m,other+[focal]); added=k1-k0
        unique=np.sum((H1-H0)*G)/np.trace(G)
        if added==0: assert abs(unique)<1e-10
        conditional.append(dict(marker=marker,focal=focal,condition_on='+'.join(other),added_df=added,estimable=added>0,unique_R2=unique if added else np.nan,note='descriptive; no conditional permutation test' if added else 'not identifiable'))
    print(marker,'completed',flush=True)
st=pd.DataFrame(rows);st['q_BH']=np.nan
for _,ix in st.groupby(['model','method','restriction']).groups.items():st.loc[ix,'q_BH']=multipletests(st.loc[ix,'p'],method='fdr_bh')[1]
st.to_csv(OUT/'model_statistics.csv',index=False)
pd.DataFrame(conditional).to_csv(OUT/'conditional_identifiability.csv',index=False)
pd.DataFrame(design).to_csv(OUT/'permutation_diagnostics.csv',index=False)
pd.DataFrame(scores).to_csv(OUT/'dbrda_scores.csv',index=False)
m_all.to_csv(OUT/'sample_manifest.csv',index=False)
pd.crosstab([m_all.marker,m_all.site,m_all.sample_date],m_all.management_group).to_csv(OUT/'design_crosstab.csv')
for label in ['management','date']:
    fig,axs=plt.subplots(1,3,figsize=(15,5),constrained_layout=True)
    for ax,marker in zip(axs,['16S','ITS','AMF']):
        m,x,y,pct,ylab=caches[(marker,label)]
        col='management_group' if label=='management' else 'sample_date'
        for group in sorted(m[col].unique()):
            ix=m[col].eq(group);ax.scatter(x[ix],y[ix],label=group,alpha=.8)
        r=st[(st.marker==marker)&(st.model==label)&(st.method=='dbRDA_Lingoes')&(st.restriction=='site_and_depth_coverage')].iloc[0]
        ax.set(title=f'{marker}: R²={r.R2:.3f}, q={r.q_BH:.4f}',xlabel=f'dbRDA 1 ({pct[0]:.1f}% of total)',ylabel=f'{ylab} ({pct[1]:.1f}% of total)')
        ax.legend(fontsize=8)
    fig.suptitle(f'{label.capitalize()}-only dbRDA; whole-profile permutations within site and depth coverage\nSeparate models retain management–date confounding',fontsize=12)
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'dbrda_{label}.{ext}',dpi=220)
    plt.close(fig)
print(st[(st.model=='management')&(st.method=='PERMANOVA_uncorrected_Bray')].to_string(index=False))
