from common import ROOT, load_counts, load_metadata, load_alpha
from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy.special import gammaln
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
O=ROOT/'reproduced/rarefied';R=O/'rarefaction';R.mkdir(parents=True,exist_ok=True)
(O/'audit').mkdir(exist_ok=True)
D=['16S','ITS','AMF'];lab={'16S':'16S','ITS':'ITS','AMF':'18S (AMF)'};B=200;SEED=20260926
colors={'NT':'#168a76','vtill':'#d76b13'}
def expected(v,k):
 v=v[v>0];N=int(v.sum())
 if k==0:return 0.
 if k==N:return float(len(v))
 valid=N-v>=k;l=np.full(len(v),-np.inf)
 l[valid]=gammaln(N-v[valid]+1)-gammaln(N-v[valid]-k+1)-gammaln(N+1)+gammaln(N-k+1)
 return float((-np.expm1(np.minimum(l,0))).sum())
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
print('Completed rarefaction; numerical outputs generated locally. No figures are created.')
