from common import ROOT, load_counts, load_metadata, load_alpha
"""Manuscript-only Lingoes Bray-Curtis dbRDA and four-set adjusted variance partitioning."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist,squareform
from scipy.linalg import svd
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
threadpool_limits(limits=1)
P=ROOT/'reproduced/full_dbrda';P.mkdir(parents=True,exist_ok=True)
B=9999;SEED=20260923
names=['Management','Date','Site','Depth'];cols=['management_group','sample_date','site','Depth_cm']
manifest=load_metadata()
models=[];fractions=[];scores=[];summary=[];unique=[];cache={}
def projection(X):
 X=X-X.mean(0);u,s,_=svd(X,full_matrices=False);q=u[:,s>1e-9];return q@q.T,q.shape[1]
def coords(d):
 A=-.5*d*d;G=A-A.mean(0)-A.mean(1)[:,None]+A.mean();c=max(0,-np.linalg.eigvalsh(G)[0]);d2=d*d+2*c;np.fill_diagonal(d2,0)
 A=-.5*d2;G=A-A.mean(0)-A.mean(1)[:,None]+A.mean();e,u=np.linalg.eigh(G);ok=e>1e-10;Y=u[:,ok]*np.sqrt(e[ok]);assert np.allclose(squareform(pdist(Y))**2,d2,atol=1e-9);return Y,c
for di,domain in enumerate(['16S','ITS','AMF']):
 m=manifest[manifest.marker==domain].reset_index(drop=True);n=len(m)
 expected=load_alpha(domain).index;assert set(expected)==set(m.sample_id)
 info,a=load_counts(domain);assert not info['taxa_are_rows'][0]
 lookup={s:i for i,s in enumerate(info['dimnames'][0])};v=a[[lookup[s] for s in m.sample_id]].astype(float);v/=v.sum(1)[:,None]
 Y,c=coords(squareform(pdist(v,metric='braycurtis')));total=np.sum(Y*Y);G=Y@Y.T
 blocks=[pd.get_dummies(m[x],drop_first=True,dtype=float).to_numpy() for x in cols];fits={0:(0.,0.,0,np.zeros((n,n)))}
 for mask in range(1,16):
  X=np.column_stack([blocks[i] for i in range(4) if mask&(1<<i)]);H,rank=projection(X);r2=np.trace(H@G)/total;adj=1-(1-r2)*(n-1)/(n-rank-1);fits[mask]=(r2,adj,rank,H)
  models.append(dict(marker=domain,mask=mask,terms='+'.join(names[i] for i in range(4) if mask&(1<<i)),df=rank,R2=r2,adjusted_R2=adj))
 r2,adj,rank,H=fits[15];fitted=H@Y;u,s,vt=svd(fitted,full_matrices=False);xy=Y@vt[:2].T;pct=s[:2]**2/total*100
 # Site scores project observed communities onto constrained axes (include residual scatter).
 assert np.isclose(np.sum(s*s),r2*total);assert np.allclose(H@H,H)
 F=(r2/rank)/((1-r2)/(n-rank-1))
 profiles=[g.index.to_numpy() for _,g in m.groupby('Site_Name')]
 strata={}
 for ids in profiles:
  signature=tuple(sorted(m.loc[ids,'Depth_cm']));strata.setdefault(signature,[]).append(ids)
 rng=np.random.default_rng(SEED+di);exceed=0
 for _ in range(B):
  ix=np.arange(n)
  for groups in strata.values():
   order=rng.permutation(len(groups))
   for target,src in zip(groups,order):ix[target]=rng.permutation(groups[src])
  assert len(set(ix))==n
  ss=np.sum(H*G[np.ix_(ix,ix)]);fp=(ss/rank)/((total-ss)/(n-rank-1));exceed+=fp>=F-1e-12
 p=(exceed+1)/(B+1)
 summary.append(dict(marker=domain,n=n,profiles=len(profiles),model_df=rank,residual_df=n-rank-1,R2=r2,adjusted_R2=adj,pseudo_F=F,p=p,permutations=B,seed=SEED+di,lingoes_constant=c,dbRDA1_total_percent=pct[0],dbRDA2_total_percent=pct[1]))
 for i,name in enumerate(names):
  other=15^(1<<i);increment=rank-fits[other][2];value=adj-fits[other][1]
  unique.append(dict(marker=domain,component=name+' unique',adjusted_fraction=value,added_df=increment,estimable=increment>0))
 # Solve all 15 exclusive Venn fractions: each submodel explains fractions intersecting its predictors.
 masks=list(range(1,16));A=np.array([[int(bool(i&j)) for j in masks] for i in masks]);z=np.linalg.solve(A,np.array([fits[i][1] for i in masks]));assert np.allclose(A@z,[fits[i][1] for i in masks]);assert np.isclose(z.sum(),adj)
 for mask,val in zip(masks,z):fractions.append(dict(marker=domain,mask=mask,component=' & '.join(names[i] for i in range(4) if mask&(1<<i)),adjusted_fraction=val))
 uq=sum(adj-fits[15^(1<<i)][1] for i in range(4));unique.extend([dict(marker=domain,component='Shared (combined)',adjusted_fraction=adj-uq,added_df=np.nan,estimable=False),dict(marker=domain,component='Unexplained',adjusted_fraction=1-adj,added_df=np.nan,estimable=False)])
 for i,row in m.iterrows():scores.append(dict(marker=domain,sample_id=row.sample_id,sample_date=row.sample_date,management=row.management_group,site=row.site,depth=row.Depth_cm,dbRDA1=xy[i,0],dbRDA2=xy[i,1]))
 cache[domain]=(m,xy,pct)
 print(domain,'R2',r2,'adj',adj,'p',p,flush=True)
st=pd.DataFrame(summary);st['q_BH']=multipletests(st.p,method='fdr_bh')[1];st.to_csv(P/'dbrda_model_statistics.csv',index=False)
for fn,data in [('variance_partition_all_fractions',fractions),('variance_partition_summary',unique),('subset_models',models),('dbrda_scores',scores)]:pd.DataFrame(data).to_csv(P/f'{fn}.csv',index=False)
manifest[['marker','sample_id','soil_uid','sample_date','management_group','site','Depth_cm']].to_csv(P/'sample_manifest.csv',index=False)
plt.rcParams.update({'font.size':11,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
dates=sorted(manifest.sample_date.unique());colors=dict(zip(dates,['#6554a4','#1b9e77','#e6ab02','#c23b32']));shapes={'NT':'o','till':'^'}
fig,axs=plt.subplots(1,3,figsize=(16,7));fig.subplots_adjust(left=.06,right=.98,top=.79,bottom=.31,wspace=.3)
for ax,domain in zip(axs,cache):
 m,xy,pct=cache[domain];r=st.set_index('marker').loc[domain]
 for mg in shapes:
  ix=m.management_group.eq(mg);ax.scatter(xy[ix,0],xy[ix,1],c=[colors[d] for d in m.loc[ix,'sample_date']],marker=shapes[mg],s=44,edgecolors='white',linewidths=.4,alpha=.9)
 ax.axhline(0,color='#ddd',lw=.7);ax.axvline(0,color='#ddd',lw=.7)
 ax.set(title=f'{domain} (n = {len(m)})',xlabel=f'dbRDA 1 ({pct[0]:.1f}% of total)',ylabel=f'dbRDA 2 ({pct[1]:.1f}% of total)');ax.set_aspect('equal',adjustable='datalim')
 ax.text(.5,-.23,f'Full model R² = {r.R2:.3f}; adjusted R² = {r.adjusted_R2:.3f}\nPseudo-F = {r.pseudo_F:.2f}; p = {r.p:.4f}; BH q = {r.q_BH:.4f}',ha='center',va='top',transform=ax.transAxes,fontsize=10,bbox=dict(facecolor='#f1f4f6',edgecolor='#ddd',boxstyle='round,pad=.5'))
fig.suptitle('dbRDA of manuscript microbial communities',fontsize=20,y=.97);fig.text(.5,.9,'Management + sampling date + site code + depth | Lingoes-corrected Bray–Curtis',ha='center',fontsize=12)
handles=[Line2D([],[],marker='o',ls='',color=colors[d],label=d) for d in dates]+[Line2D([],[],marker=s,ls='',color='#555',label=mg) for mg,s in shapes.items()]
fig.legend(handles=handles,loc='lower center',ncol=6,bbox_to_anchor=(.5,.065),frameon=False)
fig.text(.5,.035,'Exploratory full-model test: 9,999 permutations preserving profile membership; BH correction across markers.',ha='center',fontsize=10)
fig.text(.5,.008,'Management is aliased with date; full-model significance does not identify a unique tillage effect.',ha='center',fontsize=10,color='#8b4513')
for ext in ['png','pdf','svg']:fig.savefig(P/f'dbrda_manuscript.{ext}',dpi=240)
plt.close(fig)
uq=pd.DataFrame(unique);cats=['Management unique','Date unique','Site unique','Depth unique','Shared (combined)','Unexplained'];barcolors=['#bcbcbc','#6554a4','#267dab','#1b9e77','#e5ab35','#c8cdd2']
fig,axs=plt.subplots(1,3,figsize=(16,7),sharex=True);fig.subplots_adjust(left=.15,right=.94,top=.79,bottom=.23,wspace=.24)
for ax,domain in zip(axs,cache):
 vals=uq[uq.marker==domain].set_index('component').loc[cats,'adjusted_fraction'].to_numpy()*100
 ax.barh(np.arange(6),vals,color=barcolors,height=.6);ax.axvline(0,color='#555',lw=.8);ax.set_yticks(np.arange(6));ax.set_yticklabels(cats if domain=='16S' else ['']*6);ax.invert_yaxis();ax.set(title=domain,xlabel='Adjusted explained variation (%)',xlim=(-6,100))
 for j,value in enumerate(vals):
  label='Not separable from date' if j==0 else f'{value:.1f}%'
  ax.text(max(value,0)+1.4,j,label,va='center',fontsize=9.5)
fig.suptitle('Variance partitioning: what can be separated?',fontsize=20,y=.96)
fig.text(.5,.88,'Four predictor sets: management, sampling date, site code, and depth',ha='center',fontsize=12)
fig.text(.5,.13,'Unique fractions control for all other sets. Shared combines all overlaps; it cannot be assigned to one predictor.',ha='center',fontsize=10)
fig.text(.5,.085,'Management adds zero degrees of freedom after date: its unique effect is not identifiable.',ha='center',fontsize=11,color='#8b4513')
fig.text(.5,.045,'Adjusted fractions may be negative. They are retained, not interpreted as negative biological variation.',ha='center',fontsize=10)
for ext in ['png','pdf','svg']:fig.savefig(P/f'variance_partitioning_manuscript.{ext}',dpi=240)
plt.close(fig)
(P/'README.md').write_text('''# Manuscript dbRDA and variance partitioning

Cohort: exactly the manuscript MF NT/vtill cohorts, 16S n=44, ITS n=43, AMF n=44; identities verified against manuscript/audit/*_sample_alpha.csv. Uses the retained RDS count exports, explicitly reordered by sample ID. Relative abundances, no additional rarefaction or feature filtering. Dates and depth intervals are from the control workbook. Site is the original two-level manuscript site code, not each physical profile or geographic coordinates. Carbon is omitted because resolving its metadata discrepancies was outside this analysis.

Model: community ~ management + categorical sample_date + site + categorical depth. Bray-Curtis distances receive Lingoes correction. dbRDA is computed as multivariate least-squares projection of all positive corrected principal coordinates onto the centered predictor space. SVD handles rank deficiency. The fitted response defines constrained axes; plotted scores project observed communities onto these axes (not just fitted centroids). Axis percentages use total corrected inertia. Python/NumPy/SciPy implement the documented linear algebra; R/vegan was not available and these results were not cross-validated in vegan.

Full model has five effective predictor degrees of freedom, rather than six: till is identical to the November 28 date indicator. Full-model pseudo-F tests use 9,999 response permutations: profiles exchange only with profiles having the same depth coverage, while observations may also exchange within profile. This preserves paired profile membership and permits a joint test including depth. Results assume exchangeability of profiles within coverage strata and exchangeability of depths under the full null. Broader field dependence is not controlled; p-values are exploratory. These tests do not identify a unique management or date effect. BH correction is across three full-model tests. No individual predictor p-values are inferred from variance fractions.

Variance partitioning uses adjusted R² = 1 - (1-R²)*(n-1)/(n-rank-1), computed for all 15 combinations of the four sets. Unique fractions equal adjusted R²(full) minus adjusted R²(without the focal set). Shared combined equals adjusted R²(full) minus the sum of unique fractions. All 15 exclusive overlap components are also supplied; their linear sums reconstruct all subset models. Unexplained = 1 - adjusted R²(full). Negative fractions are retained. The unique management value is algebraically zero because it adds no rank after date; this is non-identifiability, not evidence of no management association. Fraction estimates are descriptive and use sample-level adjustment despite profile dependence.

Validation: sample-set equality; corrected coordinate distance reconstruction; projection idempotence; fitted inertia equals explained SS; all permutation mappings are bijections; all partition fractions reconstruct each model and sum to full adjusted R².

Outputs: PNG/PDF/SVG figures; model statistics; ordination scores; all subset model fits; all 15 partition fractions; summary fractions; sample manifest. Seeds and correction constants are recorded in dbrda_model_statistics.csv.

Reproduce from workspace root: python analysis_outputs/manuscript_dbrda/run_analysis.py

Methods: https://vegandevs.github.io/vegan/reference/dbrda.html ; https://vegandevs.github.io/vegan/reference/varpart.html ; https://vegandevs.github.io/vegan/reference/RsquareAdj.html
''')
print(st.to_string(index=False));print(uq.to_string(index=False))
