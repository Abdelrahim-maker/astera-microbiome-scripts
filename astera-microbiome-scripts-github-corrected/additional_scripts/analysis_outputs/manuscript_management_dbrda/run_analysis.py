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
P=Path('analysis_outputs/manuscript_management_dbrda');P.mkdir(exist_ok=True)
B=9999;SEED=20260923
names=['Management','Site','Depth'];cols=['management_group','site','Depth_cm']
manifest=pd.read_csv('analysis_outputs/manuscript_sample_date_ordination/manuscript_sample_manifest.csv')
book=pd.read_excel('CFARM_Co_2022_sample_control_list_2023.05.12.xlsx');book=book[book.UID_code.astype(str).str.fullmatch(r'CS\.\d+')]
manifest=manifest.merge(book[['UID_code','Site_Name']],left_on='soil_uid',right_on='UID_code',validate='many_to_one')
models=[];fractions=[];scores=[];summary=[];unique=[];cache={}
def projection(X):
 X=X-X.mean(0);u,s,_=svd(X,full_matrices=False);q=u[:,s>1e-9];return q@q.T,q.shape[1]
def coords(d):
 A=-.5*d*d;G=A-A.mean(0)-A.mean(1)[:,None]+A.mean();c=max(0,-np.linalg.eigvalsh(G)[0]);d2=d*d+2*c;np.fill_diagonal(d2,0)
 A=-.5*d2;G=A-A.mean(0)-A.mean(1)[:,None]+A.mean();e,u=np.linalg.eigh(G);ok=e>1e-10;Y=u[:,ok]*np.sqrt(e[ok]);assert np.allclose(squareform(pdist(Y))**2,d2,atol=1e-9);return Y,c
for di,domain in enumerate(['16S','ITS','AMF']):
 m=manifest[manifest.marker==domain].reset_index(drop=True);n=len(m)
 expected=pd.read_csv(f'manuscript/audit/{domain}_sample_alpha.csv').iloc[:,0];assert set(expected)==set(m.sample_id)
 info=json.loads(Path(f'manuscript/audit/{domain}_matrix_metadata.json').read_text());a=np.load(f'manuscript/audit/{domain}_counts.npz')['values'];assert not info['taxa_are_rows'][0]
 lookup={s:i for i,s in enumerate(info['dimnames'][0])};v=a[[lookup[s] for s in m.sample_id]].astype(float);v/=v.sum(1)[:,None]
 Y,c=coords(squareform(pdist(v,metric='braycurtis')));total=np.sum(Y*Y);G=Y@Y.T
 blocks=[pd.get_dummies(m[x],drop_first=True,dtype=float).to_numpy() for x in cols];fits={0:(0.,0.,0,np.zeros((n,n)))}
 for mask in range(1,8):
  X=np.column_stack([blocks[i] for i in range(3) if mask&(1<<i)]);H,rank=projection(X);r2=np.trace(H@G)/total;adj=1-(1-r2)*(n-1)/(n-rank-1);fits[mask]=(r2,adj,rank,H)
  models.append(dict(marker=domain,mask=mask,terms='+'.join(names[i] for i in range(3) if mask&(1<<i)),df=rank,R2=r2,adjusted_R2=adj))
 r2,adj,rank,H=fits[1];fitted=H@Y;u,s,vt=svd(fitted,full_matrices=False)
 assert rank==1
 residual=Y-fitted;ur,sr,vtr=svd(residual,full_matrices=False)
 x=Y@vt[0];y=ur[:,0]*sr[0]
 if x[m.management_group.eq('till')].mean()<x[m.management_group.eq('NT')].mean():x=-x
 xy=np.column_stack([x,y]);pct=np.array([s[0]**2,sr[0]**2])/total*100
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
   for target,src in zip(groups,order):
    target=sorted(target,key=lambda i:m.loc[i,'Depth_cm']);source=sorted(groups[src],key=lambda i:m.loc[i,'Depth_cm']);ix[target]=source
  assert len(set(ix))==n
  ss=np.sum(H*G[np.ix_(ix,ix)]);fp=(ss/rank)/((total-ss)/(n-rank-1));exceed+=fp>=F-1e-12
 p=(exceed+1)/(B+1)
 summary.append(dict(marker=domain,n=n,profiles=len(profiles),model_df=rank,residual_df=n-rank-1,R2=r2,adjusted_R2=adj,pseudo_F=F,p=p,permutations=B,seed=SEED+di,lingoes_constant=c,dbRDA1_total_percent=pct[0],residual_PCoA1_total_percent=pct[1]))
 full_r2,full_adj,full_rank,_=fits[7]
 for i,name in enumerate(names):
  other=7^(1<<i);increment=full_rank-fits[other][2];value=full_adj-fits[other][1]
  unique.append(dict(marker=domain,component=name+' unique',adjusted_fraction=value,added_df=increment,estimable=increment>0))
 # Solve all 15 exclusive Venn fractions: each submodel explains fractions intersecting its predictors.
 masks=list(range(1,8));A=np.array([[int(bool(i&j)) for j in masks] for i in masks]);z=np.linalg.solve(A,np.array([fits[i][1] for i in masks]));assert np.allclose(A@z,[fits[i][1] for i in masks]);assert np.isclose(z.sum(),full_adj)
 for mask,val in zip(masks,z):fractions.append(dict(marker=domain,mask=mask,component=' & '.join(names[i] for i in range(3) if mask&(1<<i)),adjusted_fraction=val))
 uq=sum(full_adj-fits[7^(1<<i)][1] for i in range(3));unique.extend([dict(marker=domain,component='Shared (combined)',adjusted_fraction=full_adj-uq,added_df=np.nan,estimable=False),dict(marker=domain,component='Unexplained',adjusted_fraction=1-full_adj,added_df=np.nan,estimable=False)])
 for i,row in m.iterrows():scores.append(dict(marker=domain,sample_id=row.sample_id,sample_date=row.sample_date,management=row.management_group,site=row.site,depth=row.Depth_cm,dbRDA1=xy[i,0],residual_PCoA1=xy[i,1]))
 cache[domain]=(m,xy,pct)
 print(domain,'R2',r2,'adj',adj,'p',p,flush=True)
st=pd.DataFrame(summary);st['q_BH']=multipletests(st.p,method='fdr_bh')[1];st.to_csv(P/'dbrda_model_statistics.csv',index=False)
for fn,data in [('variance_partition_all_fractions',fractions),('variance_partition_summary',unique),('subset_models',models),('dbrda_scores',scores)]:pd.DataFrame(data).to_csv(P/f'{fn}.csv',index=False)
manifest[['marker','sample_id','soil_uid','sample_date','management_group','site','Depth_cm']].to_csv(P/'sample_manifest.csv',index=False)
plt.rcParams.update({'font.size':11,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
dates=sorted(manifest.sample_date.unique());colors={'NT':'#247ba0','till':'#d87532'};shapes={'NT':'o','till':'^'}
fig,axs=plt.subplots(1,3,figsize=(16,7));fig.subplots_adjust(left=.06,right=.98,top=.79,bottom=.31,wspace=.3)
for ax,domain in zip(axs,cache):
 m,xy,pct=cache[domain];r=st.set_index('marker').loc[domain]
 for mg in shapes:
  ix=m.management_group.eq(mg);ax.scatter(xy[ix,0],xy[ix,1],color=colors[mg],marker=shapes[mg],s=44,edgecolors='white',linewidths=.4,alpha=.9)
 ax.axhline(0,color='#ddd',lw=.7);ax.axvline(0,color='#ddd',lw=.7)
 ax.set(title=f'{domain} (n = {len(m)})',xlabel=f'dbRDA 1 ({pct[0]:.1f}% of total)',ylabel=f'Residual PCoA 1 ({pct[1]:.1f}% of total)');ax.set_aspect('equal',adjustable='datalim')
 ax.text(.5,-.23,f'Management R² = {r.R2:.3f}; adjusted R² = {r.adjusted_R2:.3f}\nPseudo-F = {r.pseudo_F:.2f}; p = {r.p:.4f}; BH q = {r.q_BH:.4f}',ha='center',va='top',transform=ax.transAxes,fontsize=10,bbox=dict(facecolor='#f1f4f6',edgecolor='#ddd',boxstyle='round,pad=.5'))
fig.suptitle('Management-constrained dbRDA: NT versus till',fontsize=20,y=.97);fig.text(.5,.9,'Community ~ management | Lingoes-corrected Bray–Curtis | manuscript samples only',ha='center',fontsize=12)
handles=[Line2D([],[],marker=shapes[mg],ls='',color=colors[mg],label=mg) for mg in shapes]
fig.legend(handles=handles,loc='lower center',ncol=6,bbox_to_anchor=(.5,.065),frameon=False)
fig.text(.5,.035,'Exploratory management test: 9,999 permutations preserving profile membership; BH correction across markers.',ha='center',fontsize=10)
fig.text(.5,.008,'Only one constrained axis exists; the vertical axis is residual. Site, depth, and date are not adjusted here.',ha='center',fontsize=10,color='#8b4513')
for ext in ['png','pdf','svg']:fig.savefig(P/f'dbrda_manuscript.{ext}',dpi=240)
plt.close(fig)
uq=pd.DataFrame(unique);cats=['Management unique','Site unique','Depth unique','Shared (combined)','Unexplained'];barcolors=['#d87532','#267dab','#1b9e77','#e5ab35','#c8cdd2']
fig,axs=plt.subplots(1,3,figsize=(16,7),sharex=True);fig.subplots_adjust(left=.15,right=.94,top=.79,bottom=.23,wspace=.24)
for ax,domain in zip(axs,cache):
 vals=uq[uq.marker==domain].set_index('component').loc[cats,'adjusted_fraction'].to_numpy()*100
 ax.barh(np.arange(5),vals,color=barcolors,height=.6);ax.axvline(0,color='#555',lw=.8);ax.set_yticks(np.arange(5));ax.set_yticklabels(cats if domain=='16S' else ['']*5);ax.invert_yaxis();ax.set(title=domain,xlabel='Adjusted explained variation (%)',xlim=(-6,100))
 for j,value in enumerate(vals):
  label=f'{value:.1f}%'
  ax.text(max(value,0)+1.4,j,label,va='center',fontsize=9.5)
fig.suptitle('Variance partitioning: management, site, and depth',fontsize=20,y=.96)
fig.text(.5,.88,'Three predictor sets | adjusted R² | sampling date excluded',ha='center',fontsize=12)
fig.text(.5,.13,'Unique fractions control for all other sets. Shared combines all overlaps; it cannot be assigned to one predictor.',ha='center',fontsize=10)
fig.text(.5,.085,'Management unique means adjusted for site and depth only; sampling date remains confounded.',ha='center',fontsize=11,color='#8b4513')
fig.text(.5,.045,'Adjusted fractions may be negative. They are retained, not interpreted as negative biological variation.',ha='center',fontsize=10)
for ext in ['png','pdf','svg']:fig.savefig(P/f'variance_partitioning_manuscript.{ext}',dpi=240)
plt.close(fig)
(P/'README.md').write_text("""# Management-focused dbRDA and variance partitioning

Exactly the manuscript MF cohorts: 16S 44, ITS 43, AMF 44. Counts use retained manuscript/audit exports aligned by sample ID and verified against manuscript sample lists. Relative abundances and Lingoes-corrected Bray-Curtis distances are used consistently for ordination and partitioning. Date/depth/profile identifiers originate in the control workbook; site is the two-level manuscript site code.

## Management-only dbRDA
Model community ~ management. A two-level predictor provides exactly ONE constrained axis. The horizontal axis projects observed communities onto that constrained direction; the vertical axis is the first residual PCoA axis (not a second dbRDA axis). Axis percentages refer to total corrected inertia. The management effect has no adjustment for site, depth, or date. Positive horizontal scores point toward till.

Management pseudo-F uses 9,999 permutations of whole profiles, exchanged only among profiles with identical observed depth coverage, matching depths within exchanged profiles. BH correction across three marker tests. This is exploratory and assumes exchangeability of profiles within depth-coverage strata; broader field clustering is not controlled. It differs from original manuscript unrestricted 999-permutation tests, and corrected distances alter R² and F for ITS and AMF, so statistics need not match the original PERMANOVA. No statistical equivalence is inferred from nonsignificance.

## Three-set variance partitioning
Sets: management, site, depth. Date deliberately excluded. Each set is categorical. All seven subset models are fitted by orthogonal projection of full positive Lingoes-corrected principal coordinates. Adjusted R² = 1-(1-R²)*(n-1)/(n-rank-1). Unique fractions are full adjusted R² minus adjusted R² of the model without the focal set. Shared combines all exclusive overlaps. All seven individual fractions and subset fits are exported. Negative adjusted fractions are retained. Descriptive fractions have no individual significance tests. Sample-level adjusted R² does not itself account for profile dependence.

Management unique is conditional ONLY on site and depth; it cannot be interpreted as independent of sampling date. Till and the November 28 date indicator are identical. These models estimate management-associated variation, not causal tillage effects. Carbon is outside this model.

Implementation: NumPy/SciPy dbRDA linear algebra (R/vegan unavailable). Validation checks sample identities, distance reconstruction, projector idempotence, fitted sums of squares, rank-one management constraint, bijective permutations, and reconstruction of every subset model from exclusive fractions. Not cross-validated in vegan.

Reproduce: python analysis_outputs/manuscript_management_dbrda/run_analysis.py

References: https://vegandevs.github.io/vegan/reference/dbrda.html ; https://vegandevs.github.io/vegan/reference/varpart.html
""")
print(st.to_string(index=False));print(uq.to_string(index=False))
