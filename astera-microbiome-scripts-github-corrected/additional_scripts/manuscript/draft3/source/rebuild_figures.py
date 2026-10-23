from pathlib import Path
import json, shutil
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
threadpool_limits(1)
O=Path('manuscript/draft3'); A=Path('analysis_outputs'); D=['16S','ITS','AMF']; labels={'16S':'16S','ITS':'ITS','AMF':'18S (AMF)','Joint':'Cross-marker','cross':'Cross-marker'}
plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
colors=['#168a76','#d76b13']; rng=np.random.default_rng(7)
def save(fig,name):
 for ext in ['png','pdf','svg']:fig.savefig(O/'figures'/f'{name}.{ext}',dpi=300,bbox_inches='tight')
 plt.close(fig)
def read(p):return pd.read_csv(A/p)
alpha=pd.concat([read(f'{d}/{d}_alpha_wilcox.csv') for d in D],ignore_index=True);alpha['q_BH']=multipletests(alpha.p_value,method='fdr_bh')[1];alpha.to_csv(O/'audit/alpha_tests_with_BH.csv',index=False)
fig,axs=plt.subplots(2,3,figsize=(11,7))
for j,d in enumerate(D):
 a=pd.read_csv(f'manuscript/audit/{d}_sample_alpha.csv')
 for i,k in enumerate(['Observed','Shannon']):
  ax=axs[i,j];vals=[a.loc[a.management==t,k] for t in ['NT','vtill']];bp=ax.boxplot(vals,patch_artist=True,showfliers=False,medianprops={'color':'black','linewidth':1.7})
  for b,c in zip(bp['boxes'],colors):b.set(facecolor=c,alpha=.25)
  for x,(v,c) in enumerate(zip(vals,colors),1):ax.scatter(x+rng.uniform(-.12,.12,len(v)),v,s=17,color=c,alpha=.7)
  q=alpha[(alpha.domain==d)&(alpha.metric==k)].iloc[0];ax.set(title=f'{chr(65+i*3+j)}  {labels[d]}\np = {q.p_value:.3f}; BH q = {q.q_BH:.3f}',ylabel='Observed ASVs' if i==0 else 'Shannon entropy',xticks=[1,2],xticklabels=['NT','Till'])
fig.tight_layout();save(fig,'figure1_alpha')
# Reconstruct matrices and verify all alpha quantities against supplied count exports.
mat={}; meta={}; audit=[]
for d in D:
 info=json.loads(Path(f'manuscript/audit/{d}_matrix_metadata.json').read_text());a=np.load(f'manuscript/audit/{d}_counts.npz')['values'];allids=info['dimnames'][0]
 m=pd.read_csv(f'manuscript/audit/{d}_sample_alpha.csv',index_col=0);v=a[[allids.index(s) for s in m.index]];mat[d]=pd.DataFrame(v,index=m.index);meta[d]=m
 rel=v/v.sum(1)[:,None];h=-(np.where(rel>0,rel*np.log(np.maximum(rel,1e-300)),0)).sum(1)
 assert np.allclose(v.sum(1),m.reads) and np.allclose((v>0).sum(1),m.Observed) and np.allclose(h,m.Shannon)
 audit.append({'marker':d,'n':len(m),'median_retained_reads':m.reads.median(),'ASV_columns':v.shape[1],'nonzero_ASVs':int((v.sum(0)>0).sum())})
pd.DataFrame(audit).to_csv(O/'audit/count_validation.csv',index=False)
# Three-marker PCoA plus cross-marker standardized CLR-PCA.
co=read('manuscript_sample_date_ordination/ordination_coordinates.csv');summ=read('manuscript_sample_date_ordination/ordination_summary.csv').set_index('marker')
fig,axs=plt.subplots(2,2,figsize=(10,9))
for ax,d,letter in zip(axs.flat,D,'ABC'):
 s=co[co.marker==d];r=summ.loc[d];test=read(f'{d}/{d}_nmds_permanova.csv').query("term == 'Model'").iloc[0]
 for mg,c,shape in zip(['NT','till'],colors,['o','^']):
  x=s[s.management_group==mg];ax.scatter(x.PCoA1,x.PCoA2,c=c,marker=shape,s=38,alpha=.8,label=mg.upper() if mg=='NT' else 'Till')
 ax.set(title=f'{letter}  {labels[d]} (n = {len(s)})\nManagement R² = {test.R2:.3f}; p = {test["Pr(>F)"]:.3f}',xlabel=f'PCoA 1 ({r.axis1_percent:.1f}%)',ylabel=f'PCoA 2 ({r.axis2_percent:.1f}%)')
 ax.axhline(0,color='#ddd',lw=.6);ax.axvline(0,color='#ddd',lw=.6)
common=[s for s in mat['16S'].index if all(s in mat[d].index for d in D)];blocks=[]
for d in D:
 v=mat[d].loc[common].to_numpy();v=v[:,(v>0).mean(0)>=.05]
 if v.shape[1]>200:v=v[:,np.argsort(-v.sum(0),kind='stable')[:200]]
 v=np.log(v+1);v-=v.mean(1)[:,None];blocks.append(v)
X=np.column_stack(blocks);X=(X-X.mean(0))/X.std(0,ddof=1);u,s,vt=svd(X,full_matrices=False);xy=u[:,:2]*s[:2];pct=s[:2]**2/(s*s).sum()*100
mg=meta['16S'].loc[common].management;ax=axs[1,1]
for group,c,shape in zip(['NT','vtill'],colors,['o','^']):
 ix=mg.eq(group);ax.scatter(xy[ix,0],xy[ix,1],color=c,marker=shape,s=38,alpha=.8)
G=X@X.T;Z=pd.get_dummies(mg,drop_first=True,dtype=float).to_numpy();Z=Z-Z.mean(0);H=Z@np.linalg.pinv(Z);r2=np.trace(H@G)/np.trace(G)
st=read('joint_rpca_proxy_permanova.csv').query("term == 'Model'").iloc[0];assert np.isclose(r2,st.R2)
ax.set(title=f'D  Cross-marker (n = {len(common)})\nManagement R² = {st.R2:.3f}; p = {st["Pr(>F)"]:.3f}',xlabel=f'PC 1 ({pct[0]:.1f}%)',ylabel=f'PC 2 ({pct[1]:.1f}%)');ax.axhline(0,color='#ddd',lw=.6);ax.axvline(0,color='#ddd',lw=.6)
axs[0,0].legend(frameon=False);fig.tight_layout();save(fig,'figure2_composition')
pd.DataFrame({'sample_id':common,'PC1':xy[:,0],'PC2':xy[:,1]}).to_csv(O/'audit/cross_marker_pca_scores.csv',index=False)
# Date-inclusive dbRDA; display full-model rather than individual term tests.
st=read('manuscript_dbrda/dbrda_model_statistics.csv').set_index('marker');sc=read('manuscript_dbrda/dbrda_scores.csv');dates=sorted(sc.sample_date.unique());dc=dict(zip(dates,['#6554a4','#168a76','#c23b32','#d4a017']))
fig,axs=plt.subplots(1,3,figsize=(11,5))
for ax,d,letter in zip(axs,D,'ABC'):
 t=sc[sc.marker==d];r=st.loc[d]
 for mg,shape in [('NT','o'),('till','^')]:
  x=t[t.management==mg];ax.scatter(x.dbRDA1,x.dbRDA2,c=[dc[z] for z in x.sample_date],marker=shape,s=40,edgecolors='white',lw=.3)
 ax.set(title=f'{letter}  {labels[d]}\nAdj. R² = {r.adjusted_R2:.3f}\nFull-model q < 0.001',xlabel=f'dbRDA 1 ({r.dbRDA1_total_percent:.1f}% of total)',ylabel=f'dbRDA 2 ({r.dbRDA2_total_percent:.1f}% of total)');ax.axhline(0,color='#ddd',lw=.6);ax.axvline(0,color='#ddd',lw=.6)
handles=[Line2D([],[],marker='o',ls='',color=dc[d],label=d) for d in dates]+[Line2D([],[],marker=s,ls='',color='#555',label=t) for t,s in [('NT','o'),('Till','^')]]
fig.legend(handles=handles,loc='lower center',ncol=6,frameon=False);fig.tight_layout(rect=(0,.12,1,1));save(fig,'figure3_dbrda')
vp=read('manuscript_dbrda/variance_partition_summary.csv');cats=['Date unique','Site unique','Depth unique','Shared (combined)']
fig,axs=plt.subplots(1,3,figsize=(11,5),sharex=True)
for ax,d,letter in zip(axs,D,'ABC'):
 v=vp[vp.marker==d].set_index('component').adjusted_fraction;vals=v.loc[cats]*100
 ax.barh(range(4),vals,color=['#6554a4','#267dab','#168a76','#d4a017'],height=.6);ax.set(yticks=range(4),yticklabels=cats if d=='16S' else ['']*4,xlabel='Adjusted fraction (%)',title=f'{letter}  {labels[d]}\nUnexplained: {v.Unexplained*100:.1f}%',xlim=(0,8.5));ax.invert_yaxis()
 for i,x in enumerate(vals):ax.text(x+.13,i,f'{x:.1f}',va='center')
fig.text(.5,.02,'Unique management is not identifiable after accounting for date.',ha='center',fontsize=12);fig.tight_layout(rect=(0,.08,1,1));save(fig,'figure4_partition')
fig,ax=plt.subplots(figsize=(7,4));x=np.arange(3)
for offset,cut,color in [(-.16,.05,'#334e68'),(.16,.25,'#a6bbc8')]:
 vals=[int((read(f'{d}/{d}_maaslin2_nt_vs_till_results.csv').qval<=cut).sum()) for d in D];bars=ax.bar(x+offset,vals,width=.28,color=color,label=f'q ≤ {cut:.2f} (total)');ax.bar_label(bars,padding=3)
ax.set(xticks=x,xticklabels=[labels[d] for d in D],ylabel='Treatment-associated ASVs',ylim=(0,53));ax.legend(frameon=False);fig.tight_layout();save(fig,'figure5_features')
net=read('network_stats_summary.csv');fig,axs=plt.subplots(1,3,figsize=(11,4))
for ax,k,label in zip(axs,['nodes','mean_degree','density'],['Connected ASVs','Mean degree','Density']):
 for i,(t,c) in enumerate(zip(['NT','till'],colors)):
  vals=[net.loc[net.label==f'{d} {t}',k].iloc[0] for d in D+['cross']];ax.bar(np.arange(4)+(i-.5)*.34,vals,width=.3,color=c,label='Till' if t=='till' else t)
 ax.set(xticks=range(4),xticklabels=[labels[d] for d in D+['cross']],ylabel=label);ax.tick_params(axis='x',labelrotation=20)
axs[0].legend(frameon=False);fig.tight_layout();save(fig,'figure6_networks')
b=read('imbalance_sensitivity/ordination_balanced_subsamples_permanova_permdisp.csv');j=read('imbalance_sensitivity/joint_balanced_subsamples_permanova.csv');fig,axs=plt.subplots(1,2,figsize=(10,4))
for ax,k,jk,label in zip(axs,['permanova_r2','permanova_p'],['r2','p_value'],['PERMANOVA R²','Permutation p value']):
 vals=[b.loc[b.domain==d,k].values for d in D]+[j[jk].values];ax.boxplot(vals,showfliers=False,medianprops={'color':'black','linewidth':1.8})
 for i,v in enumerate(vals,1):ax.scatter(i+rng.uniform(-.1,.1,len(v)),v,s=20,color='#334e68',alpha=.7)
 ax.set(xticks=range(1,5),xticklabels=[labels[d] for d in D+['Joint']],ylabel=label);ax.tick_params(axis='x',labelrotation=15)
axs[1].axhline(.05,color='#b32632',ls='--');fig.tight_layout();save(fig,'figure7_sensitivity')
c=read('carbon_analysis/carbon_permanova_by_domain.csv');j=read('carbon_analysis/carbon_joint_permanova.csv');j['domain']='Joint';c=pd.concat([c,j]);c=c[c.test_type.eq('sequential') & c.term.isin(['trt2','C','trt2:C'])];c.to_csv(O/'audit/carbon_sequential.csv',index=False)
fig,ax=plt.subplots(figsize=(7,4));vals=np.array([[c.loc[(c.domain==d)&(c.term==t),'R2'].iloc[0]*100 for t in ['trt2','C','trt2:C']] for d in D+['Joint']]);im=ax.imshow(vals,cmap='Blues',vmin=0,vmax=7,aspect='auto')
for i,d in enumerate(D+['Joint']):
 for k,t in enumerate(['trt2','C','trt2:C']):
  p=c.loc[(c.domain==d)&(c.term==t),'Pr(>F)'].iloc[0];ax.text(k,i,f'{vals[i,k]:.1f}%\np = {p:.3f}',ha='center',va='center',color='white' if vals[i,k]>4 else 'black')
ax.set(xticks=range(3),xticklabels=['Management','C after\nmanagement','Management\n× C'],yticks=range(4),yticklabels=[labels[d] for d in D+['Joint']]);fig.colorbar(im,ax=ax,label='Sequential explained variation (%)');fig.tight_layout();save(fig,'figure8_carbon')
# Supplemental date view with the same coordinates as Figure 2.
fig,axs=plt.subplots(1,3,figsize=(11,5))
for ax,d in zip(axs,D):
 t=co[co.marker==d];r=summ.loc[d]
 for mg,shape in [('NT','o'),('till','^')]:
  x=t[t.management_group==mg];ax.scatter(x.PCoA1,x.PCoA2,c=[dc[z] for z in x.sample_date],marker=shape,s=40)
 ax.set(title=labels[d],xlabel=f'PCoA 1 ({r.axis1_percent:.1f}%)',ylabel=f'PCoA 2 ({r.axis2_percent:.1f}%)')
fig.legend(handles=handles,loc='lower center',ncol=6,frameon=False);fig.tight_layout(rect=(0,.12,1,1));save(fig,'figureS54_date')
print('Figures rebuilt; count-derived alpha quantities and cross-marker R² verified.')
