from common import ROOT, load_counts, load_metadata, load_alpha
from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy.special import gammaln
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
O=ROOT/'reproduced/rarefied';R=O/'rarefaction';R.mkdir(parents=True,exist_ok=True)
(O/'figures').mkdir(exist_ok=True);(O/'audit').mkdir(exist_ok=True)
D=['16S','ITS','AMF'];lab={'16S':'16S','ITS':'ITS','AMF':'18S (AMF)'};B=200;SEED=20260926
colors={'NT':'#168a76','vtill':'#d76b13'}
def expected(v,k):
 v=v[v>0];N=int(v.sum())
 if k==0:return 0.
 if k==N:return float(len(v))
 valid=N-v>=k;l=np.full(len(v),-np.inf)
 l[valid]=gammaln(N-v[valid]+1)-gammaln(N-v[valid]-k+1)-gammaln(N+1)+gammaln(N-k+1)
 return float((-np.expm1(np.minimum(l,0))).sum())
def save(fig,name):
 for ext in ['png','pdf','svg']:fig.savefig(O/'figures'/f'{name}.{ext}',dpi=300,bbox_inches='tight')
 plt.close(fig)
# Small exact enumeration check of analytic expectation.
assert np.isclose(expected(np.array([2,1]),2),5/3)
rows=[];depths=[];curve=[];exclusions=[];testrows=[]
for di,d in enumerate(D):
 info,raw=load_counts(d);ids=info['dimnames'][0]
 m=load_alpha(d).sort_index();v=raw[[ids.index(s) for s in m.index]]
 assert (v>=0).all() and np.array_equal(v,np.floor(v));v=v.astype(np.int64);total=v.sum(1)
 primary=int(total.min());high=int(np.sort(total)[int(np.floor(.05*len(total)))])
 for si,(scenario,depth) in enumerate([('primary_all_samples',primary),('sensitivity_95pct',high)]):
  keep=total>=depth;rng=np.random.default_rng(SEED+100*di+si)
  depths.append(dict(marker=d,scenario=scenario,depth=depth,n_retained=int(keep.sum()),n_excluded=int((~keep).sum()),n_NT=int(((m.management=='NT')&keep).sum()),n_till=int(((m.management=='vtill')&keep).sum()),repeats=B,seed=SEED+100*di+si))
  for i,sid in enumerate(m.index):
   if not keep[i]:exclusions.append(dict(marker=d,scenario=scenario,sample_id=sid,reads=int(total[i]),depth=depth,management=m.iloc[i].management));continue
   x=v[i][v[i]>0];z=rng.multivariate_hypergeometric(x,depth,size=B);assert (z.sum(1)==depth).all() and (z<=x).all() and (z>=0).all()
   rich=(z>0).sum(1);prob=z/depth;shan=-(prob*np.log(np.maximum(prob,1e-300))).sum(1);ex=expected(x,depth)
   assert abs(rich.mean()-ex)<6*rich.std(ddof=1)/np.sqrt(B)+.05
   for rep in range(B):rows.append(dict(marker=d,scenario=scenario,sample_id=sid,management=m.iloc[i].management,depth=depth,replicate=rep+1,Observed=int(rich[rep]),Shannon=float(shan[rep]),expected_richness=ex))
  print(d,scenario,depth,'retained',keep.sum(),flush=True)
 for i,sid in enumerate(m.index):
  x=v[i][v[i]>0];N=int(total[i]);grid=np.unique(np.r_[0,1,np.geomspace(1,N,70).astype(int),np.linspace(1,primary,45).astype(int),primary,N]);last=-1
  for k in grid:
   e=expected(x,int(k));assert e>=last-1e-6 and e<=len(x)+1e-7;last=e
   curve.append(dict(marker=d,sample_id=sid,management=m.iloc[i].management,reads=int(k),expected_richness=e,library_total=N))
a=pd.DataFrame(rows);a.to_csv(R/'rarefaction_replicates.csv',index=False)
g=a.groupby(['marker','scenario','sample_id','management','depth'],sort=False)
s=g.agg(Observed=('Observed','mean'),Shannon=('Shannon','mean'),richness_sd=('Observed','std'),shannon_sd=('Shannon','std'),expected_richness=('expected_richness','first')).reset_index();s.to_csv(R/'sample_alpha_rarefied.csv',index=False)
for scenario in s.scenario.unique():
 part=[]
 for d in D:
  x=s[(s.marker==d)&(s.scenario==scenario)]
  for metric in ['Observed','Shannon']:
   nt=x.loc[x.management=='NT',metric];ti=x.loc[x.management=='vtill',metric];u,p=mannwhitneyu(nt,ti,alternative='two-sided',method='asymptotic',use_continuity=True)
   part.append(dict(domain=d,scenario=scenario,metric=metric,depth=int(x.depth.iloc[0]),n_nt=len(nt),n_till=len(ti),nt_median=nt.median(),till_median=ti.median(),U=u,p_value=p))
 q=multipletests([r['p_value'] for r in part],method='fdr_bh')[1]
 for r,qq in zip(part,q):r['q_BH']=qq
 testrows+=part
st=pd.DataFrame(testrows);st.to_csv(R/'alpha_tests_rarefied.csv',index=False)
primary=st[st.scenario=='primary_all_samples'];primary.to_csv(O/'audit/alpha_tests_with_BH.csv',index=False)
pd.DataFrame(depths).to_csv(R/'depth_selection.csv',index=False);pd.DataFrame(exclusions).to_csv(R/'excluded_samples.csv',index=False);c=pd.DataFrame(curve);c.to_csv(R/'rarefaction_curves.csv',index=False)
# Repeated-draw sensitivity: tests within each draw, never pooling draw observations.
rt=[]
for (scenario,rep),x in a.groupby(['scenario','replicate']):
 part=[]
 for d in D:
  y=x[x.marker==d]
  for metric in ['Observed','Shannon']:
   _,p=mannwhitneyu(y.loc[y.management=='NT',metric],y.loc[y.management=='vtill',metric],alternative='two-sided',method='asymptotic',use_continuity=True)
   part.append(dict(scenario=scenario,replicate=rep,marker=d,metric=metric,p_value=p))
 for r,q in zip(part,multipletests([r['p_value'] for r in part],method='fdr_bh')[1]):r['q_BH']=q
 rt+=part
rt=pd.DataFrame(rt);rt.to_csv(R/'per_draw_tests.csv',index=False)
summary=rt.groupby(['scenario','marker','metric']).agg(p_min=('p_value','min'),p_median=('p_value','median'),p_max=('p_value','max'),q_min=('q_BH','min'),q_max=('q_BH','max'),draws_p_lt_005=('p_value',lambda x:int((x<.05).sum())),draws_q_lt_005=('q_BH',lambda x:int((x<.05).sum()))).reset_index();summary.to_csv(R/'per_draw_test_summary.csv',index=False)
plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,axs=plt.subplots(2,3,figsize=(11,7));rng=np.random.default_rng(7)
for j,d in enumerate(D):
 x=s[(s.marker==d)&(s.scenario=='primary_all_samples')];dep=int(x.depth.iloc[0])
 for i,metric in enumerate(['Observed','Shannon']):
  ax=axs[i,j];vals=[x.loc[x.management==mg,metric] for mg in ['NT','vtill']];bp=ax.boxplot(vals,patch_artist=True,showfliers=False,medianprops={'color':'black','linewidth':1.7})
  for b,col in zip(bp['boxes'],colors.values()):b.set(facecolor=col,alpha=.25)
  for k,(v0,col) in enumerate(zip(vals,colors.values()),1):ax.scatter(k+rng.uniform(-.12,.12,len(v0)),v0,color=col,s=18,alpha=.75)
  r=primary[(primary.domain==d)&(primary.metric==metric)].iloc[0];ax.set(title=f'{chr(65+3*i+j)}  {lab[d]} · {dep:,} reads\np = {r.p_value:.3f}; q = {r.q_BH:.3f}',xticks=[1,2],xticklabels=['NT','Till'],ylabel='Mean rarefied richness' if i==0 else 'Mean rarefied Shannon')
fig.tight_layout();save(fig,'figure1_alpha')
fig,axs=plt.subplots(2,3,figsize=(11,7))
for j,d in enumerate(D):
 dep=int(primary[primary.domain==d].depth.iloc[0])
 for i in [0,1]:
  ax=axs[i,j]
  for sid,x in c[c.marker==d].groupby('sample_id'):
   x=x[(x.reads>0)] if i==0 else x[x.reads<=dep];ax.plot(x.reads.to_numpy(),x.expected_richness.to_numpy(),color=colors[x.management.iloc[0]],lw=.8,alpha=.5)
  ax.axvline(dep,color='black',ls='--',lw=1)
  if i==0:ax.set_xscale('log')
  else:ax.set_xlim(0,dep*1.04)
  ax.set(title=f'{chr(65+3*i+j)}  {lab[d]}' + (' · full libraries' if i==0 else f' · to {dep:,} reads'),xlabel='Reads sampled (log scale)' if i==0 else 'Reads sampled',ylabel='Expected ASV richness')
fig.legend(handles=[Line2D([],[],color=col,label='Till' if mg=='vtill' else mg) for mg,col in colors.items()]+[Line2D([],[],color='black',ls='--',label='Primary rarefaction depth')],loc='lower center',ncol=3,frameon=False)
fig.tight_layout(rect=(0,.07,1,1));save(fig,'figure9_rarefaction')
(R/'README.md').write_text('''# Rarefied alpha diversity
Primary depths are the minimum retained library size for each marker; no samples are excluded. A secondary, label-blind depth is the largest observed depth retaining at least 95% of each marker cohort. Each sample is subsampled without replacement 200 times with NumPy multivariate_hypergeometric. Main tests use one mean richness and one mean natural-log Shannon value per biological observation, not 200 pseudo-replicates. Two-sided asymptotic Mann–Whitney/Wilcoxon rank-sum tests use tie and continuity corrections; BH is across six tests separately for each depth scenario. Per-draw tests describe rarefaction randomness, not independent experiments. Confounding and profile dependence remain; these are exploratory unadjusted management comparisons.

Rarefaction curves use exact hypergeometric expected richness: sum_i [1 - choose(N-n_i,k)/choose(N,k)]. Zero-only features do not affect curves. No extrapolation, additional ASV filtering, or count rounding was applied. Existing filtering/singleton limitations in the supplied processed matrices remain. Read rarefaction cannot recover previously removed taxa. Seeds, excluded IDs, all 200-draw metrics, and numerical validation are exported.

References: https://vegandevs.github.io/vegan/reference/rarefy.html ; https://numpy.org/doc/stable/reference/random/generated/numpy.random.Generator.multivariate_hypergeometric.html
''')
print(st.to_string(index=False));print(summary.to_string(index=False))
