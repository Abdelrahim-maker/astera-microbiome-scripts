"""Date-only dbRDA of NT manuscript samples; run from workspace root."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
threadpool_limits(limits=1)
OUT=Path('analysis_outputs/nt_sample_date_dbrda'); OUT.mkdir(exist_ok=True)
B=9999; SEED=20260929
manifest=pd.read_csv('analysis_outputs/nt_sample_date/sample_manifest.csv')
assert manifest.management_group.eq('NT').all()
assert not manifest.duplicated(['scope','marker','sample_id']).any()
rows=[]; scores=[]; cache={}
def projector(g):
 Z=np.eye(len(np.unique(g)))[g]
 return Z@np.linalg.inv(Z.T@Z)@Z.T-np.ones((len(g),len(g)))/len(g)
def fit(G,H,rank):
 total=np.trace(G); ss=np.sum(H*G.T); residual=total-ss
 return (ss/rank)/(residual/(len(G)-rank-1)),ss/total
for scope in ['NT_all','NT_Co_NT_only']:
 for k,marker in enumerate(['16S','ITS','AMF']):
  m=manifest[(manifest.scope==scope)&(manifest.marker==marker)].reset_index(drop=True)
  info=json.loads(Path(f'manuscript/audit/{marker}_matrix_metadata.json').read_text())
  assert not info['taxa_are_rows'][0]
  a=np.load(f'manuscript/audit/{marker}_counts.npz')['values']; lookup={s:i for i,s in enumerate(info['dimnames'][0])}
  v=a[[lookup[s] for s in m.sample_id]].astype(float)
  assert np.isfinite(v).all() and (v>=0).all() and (v.sum(1)>0).all()
  v/=v.sum(1)[:,None]; d2=squareform(pdist(v,'braycurtis'))**2
  A=-.5*d2; rawG=A-A.mean(0)-A.mean(1)[:,None]+A.mean()
  c=max(0.,-np.linalg.eigvalsh(rawG)[0]); dc2=d2+2*c; np.fill_diagonal(dc2,0)
  A=-.5*dc2; G=A-A.mean(0)-A.mean(1)[:,None]+A.mean()
  e,U=np.linalg.eigh(G); Y=U[:,e>1e-10]*np.sqrt(e[e>1e-10])
  assert np.allclose(squareform(pdist(Y))**2,dc2,atol=1e-9)
  g=pd.Categorical(m.sample_date,categories=sorted(m.sample_date.unique())).codes
  rank=len(np.unique(g))-1; H=projector(g); F,R2=fit(G,H,rank)
  assert np.allclose(H@H,H)
  fitted=H@Y; u,s,vt=np.linalg.svd(fitted,full_matrices=False)
  assert np.isclose((s*s).sum(),R2*np.trace(G))
  # Independent group sum-of-squares check.
  within=sum(dc2[np.ix_(g==z,g==z)].sum()/(2*(g==z).sum()) for z in np.unique(g))
  assert np.isclose(F,((np.trace(G)-within)/rank)/(within/(len(g)-rank-1)))
  xy=Y@vt[:2].T; pct=s[:2]**2/np.trace(G)*100
  profiles=sorted(m.Site_Name.unique()); pindex={p:i for i,p in enumerate(profiles)}
  ri=np.array([pindex[p] for p in m.Site_Name]); pg=[]; signatures=[]
  for p in profiles:
   z=m[m.Site_Name==p]; assert z.sample_date.nunique()==1 and not z.Depth_cm.duplicated().any()
   pg.append(g[z.index[0]]); signatures.append('|'.join(sorted(z.Depth_cm)))
  pg=np.array(pg); signatures=np.array(signatures)
  strata=[np.where(signatures==z)[0] for z in np.unique(signatures)]
  rng=np.random.default_rng(SEED+k); exceed=0
  for _ in range(B):
   shuffled=pg.copy()
   for ids in strata: shuffled[ids]=rng.permutation(pg[ids])
   perm=shuffled[ri]; assert np.array_equal(np.bincount(g),np.bincount(perm))
   exceed+=fit(G,projector(perm),rank)[0]>=F-1e-12
  row=dict(scope=scope,marker=marker,n=len(m),profiles=len(profiles),model='community ~ categorical sampling_date',model_df=rank,residual_df=len(m)-rank-1,R2=R2,adjusted_R2=1-(1-R2)*(len(m)-1)/(len(m)-rank-1),pseudo_F=F,p=(exceed+1)/(B+1),permutations=B,seed=SEED+k,lingoes_constant=c,dbRDA1_total_percent=pct[0],dbRDA2_total_percent=pct[1])
  rows.append(row)
  z=m.copy();z['dbRDA1']=xy[:,0];z['dbRDA2']=xy[:,1];scores.append(z)
  cache[(scope,marker)]=(m,xy,pct)
  print(row,flush=True)
st=pd.DataFrame(rows)
for scope,ix in st.groupby('scope').groups.items(): st.loc[ix,'q_BH']=multipletests(st.loc[ix,'p'],method='fdr_bh')[1]
st.to_csv(OUT/'dbrda_statistics.csv',index=False);pd.concat(scores).to_csv(OUT/'dbrda_scores.csv',index=False)
manifest.to_csv(OUT/'sample_manifest.csv',index=False)
plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
dates=sorted(manifest.sample_date.unique());colors=dict(zip(dates,['#6554a4','#1b9e77','#c23b32']))
for scope in ['NT_all','NT_Co_NT_only']:
 fig,axs=plt.subplots(1,3,figsize=(14,5.5));fig.subplots_adjust(left=.06,right=.98,top=.80,bottom=.30,wspace=.30)
 for ax,marker in zip(axs,['16S','ITS','AMF']):
  m,xy,pct=cache[(scope,marker)];r=st[(st.scope==scope)&(st.marker==marker)].iloc[0]
  for date in dates:
   ix=m.sample_date.eq(date);ax.scatter(xy[ix,0],xy[ix,1],c=colors[date],s=40,edgecolor='white',linewidth=.4)
  ax.set(title=f'{marker} (n = {len(m)})',xlabel=f'dbRDA1 ({pct[0]:.1f}% of total)',ylabel=f'dbRDA2 ({pct[1]:.1f}% of total)')
  ax.axhline(0,color='#ddd',lw=.6);ax.axvline(0,color='#ddd',lw=.6)
  ax.text(.5,-.25,f'R² = {r.R2:.3f}; adjusted R² = {r.adjusted_R2:.3f}\np = {r.p:.4f}; BH q = {r.q_BH:.4f}',transform=ax.transAxes,ha='center',va='top')
 fig.suptitle('Non-tilled samples: dbRDA constrained by sampling date',fontsize=17)
 fig.text(.5,.88,'Lingoes-corrected Bray–Curtis | '+('All NT samples' if scope=='NT_all' else 'Co_NT site only'),ha='center')
 fig.legend(handles=[Line2D([],[],marker='o',ls='',color=colors[d],label=d) for d in dates],loc='lower center',bbox_to_anchor=(.5,.035),ncol=3,frameon=False)
 fig.text(.5,.01,'9,999 profile-level permutations within depth-coverage strata; exploratory date association.',ha='center',fontsize=9)
 for ext in ['png','pdf']:fig.savefig(OUT/f'dbrda_{scope}.{ext}',dpi=220)
 plt.close(fig)
print(st.to_string(index=False))
