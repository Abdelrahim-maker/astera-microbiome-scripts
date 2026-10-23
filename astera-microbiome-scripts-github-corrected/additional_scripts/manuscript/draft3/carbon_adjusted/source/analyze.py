from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist,squareform
from scipy.linalg import orth
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits
threadpool_limits(1)
O=Path('manuscript/draft3/carbon_adjusted'); B=9999
m0=pd.read_csv('analysis_outputs/manuscript_sample_date_ordination/manuscript_sample_manifest.csv')
c=pd.read_csv('manuscript/draft3/sample_tables/sample_details.csv')
b=pd.read_excel('CFARM_Co_2022_sample_control_list_2023.05.12.xlsx');b=b[b.UID_code.astype(str).str.fullmatch(r'CS\.\d+')]
m0=m0.merge(b[['UID_code','Site_Name']],left_on='soil_uid',right_on='UID_code',validate='many_to_one').merge(c[['soil_uid','carbon_conflict','C_16S','C_ITS','C_AMF']],on='soil_uid',validate='many_to_one')
community=[];reg=[];summary=[];checks=[];man=[]
def design(m,cols):
 return np.column_stack([np.ones(len(m))]+[pd.get_dummies(m[x],drop_first=True,dtype=float).to_numpy() if x!='carbon' else m[[x]].to_numpy() for x in cols])
def H(X):
 q=orth(X);return q@q.T,q.shape[1]
def coord(v):
 d=squareform(pdist(v,'braycurtis'));J=np.eye(len(d))-np.ones_like(d)/len(d);g=-.5*J@(d*d)@J
 correction=max(0,-np.linalg.eigvalsh(g)[0]);d2=d*d+2*correction;np.fill_diagonal(d2,0);g=-.5*J@d2@J
 e,u=np.linalg.eigh(g);pos=e>1e-10;y=u[:,pos]*np.sqrt(e[pos]);assert np.allclose(squareform(pdist(y))**2,d2,atol=1e-8)
 return y,correction
for marker in ['16S','ITS','AMF']:
 base=m0[m0.marker.eq(marker)].copy();base['carbon']=base['C_'+marker]
 for variant in ['recorded','exclude_conflicting']:
  m=base[base.carbon.notna() & ((~base.carbon_conflict) if variant=='exclude_conflicting' else True)].reset_index(drop=True)
  if variant=='recorded':man.append(m)
  summary.append(dict(marker=marker,variant=variant,n=len(m),profiles=m.Site_Name.nunique(),NT_n=m.management_group.eq('NT').sum(),till_n=m.management_group.eq('till').sum(),NT_median=m.loc[m.management_group.eq('NT'),'carbon'].median(),till_median=m.loc[m.management_group.eq('till'),'carbon'].median(),surface_median=m.loc[m.Depth_cm.eq('0-10'),'carbon'].median(),deep_median=m.loc[m.Depth_cm.eq('10-30'),'carbon'].median()))
  for label,formula in [('management_depth_site','carbon ~ C(management_group, Treatment(reference="NT")) + C(Depth_cm) + C(site)'),('date_depth_site','carbon ~ C(sample_date) + C(Depth_cm) + C(site)')]:
   fit=smf.ols(formula,m).fit(cov_type='cluster',cov_kwds={'groups':m.Site_Name,'use_correction':True},use_t=True)
   for term in fit.params.index:
    ci=fit.conf_int().loc[term]
    reg.append(dict(marker=marker,variant=variant,model=label,term=term,estimate=fit.params[term],se=fit.bse[term],p=fit.pvalues[term],CI_low=ci.iloc[0],CI_high=ci.iloc[1],n=len(m),clusters=m.Site_Name.nunique(),R2=fit.rsquared))
  info=json.loads(Path(f'manuscript/audit/{marker}_matrix_metadata.json').read_text());a=np.load(f'manuscript/audit/{marker}_counts.npz')['values'];look={s:i for i,s in enumerate(info['dimnames'][0])};v=a[[look[s] for s in m.sample_id]].astype(float);v/=v.sum(1)[:,None];y,lc=coord(v);total=np.sum(y*y)
  # Descriptive full cohort increments; management/date are aliased.
  terms=['management_group','sample_date','site','Depth_cm','carbon']; hf,rf=H(design(m,terms));ssef=np.sum((y-hf@y)**2)
  for term in ['carbon','management_group','Depth_cm']:
   hr,rr=H(design(m,[t for t in terms if t!=term]));sser=np.sum((y-hr@y)**2);df=rf-rr
   community.append(dict(marker=marker,variant=variant,model='date_site_depth_management_carbon_descriptive',term=term,n=len(m),profiles=m.Site_Name.nunique(),added_df=df,residual_df=len(m)-rf,unique_R2=(sser-ssef)/total,partial_R2=(sser-ssef)/sser,unique_adjusted_R2=(1-ssef/total*(len(m)-1)/(len(m)-rf))-(1-sser/total*(len(m)-1)/(len(m)-rr)),pseudo_F=((sser-ssef)/df)/(ssef/(len(m)-rf)) if df else np.nan,p=np.nan,permutations=0,lingoes_constant=lc))
  # Paired complete profiles only: profile FE removes all between-profile confounding.
  ids=m.groupby('Site_Name').filter(lambda g:len(g)==2 and set(g.Depth_cm)=={'0-10','10-30'}).index.to_numpy();mp=m.loc[ids].reset_index(drop=True);vp=v[ids];yp,lc=coord(vp);n=len(mp);tot=np.sum(yp*yp)
  mp[['carbon','Site_Name','Depth_cm']].to_csv(O/f'{marker}_{variant}_paired_metadata.csv',index=False);pd.DataFrame(yp).to_csv(O/f'{marker}_{variant}_coordinates.csv',index=False)
  xr=design(mp,['Site_Name','Depth_cm']);xf=design(mp,['Site_Name','Depth_cm','carbon']);hr,rr=H(xr);hf,rf=H(xf);df=rf-rr;assert df==1
  res=yp-hr@yp; fitted=hr@yp;ssr=np.sum(res*res);ssf=np.sum((yp-hf@yp)**2);F=((ssr-ssf)/df)/(ssf/(n-rf))
  rng=np.random.default_rng(20261003);groups=[g.index.to_numpy() for _,g in mp.groupby('Site_Name')];exceed=0
  for _ in range(B):
   ix=np.arange(n)
   for g in groups:
    if rng.integers(2):ix[g]=g[::-1]
   ys=fitted+res[ix];er=np.sum((ys-hr@ys)**2);ef=np.sum((ys-hf@ys)**2);fp=((er-ef)/df)/(ef/(n-rf));exceed+=fp>=F-1e-12
  p=(exceed+1)/(B+1)
  # Independent least-squares residual calculation.
  direct=np.sum((yp-xf@np.linalg.lstsq(xf,yp,rcond=None)[0])**2);assert np.isclose(direct,ssf)
  assert np.allclose(hf@hf,hf) and np.isclose((ssr-ssf)/ssr,1-ssf/ssr)
  checks.append(dict(marker=marker,variant=variant,check='distance reconstruction, projection, least squares and partial R2',passed=True))
  community.append(dict(marker=marker,variant=variant,model='within_profile_carbon_adjusted_depth',term='carbon',n=n,profiles=len(groups),added_df=df,residual_df=n-rf,unique_R2=(ssr-ssf)/tot,partial_R2=(ssr-ssf)/ssr,unique_adjusted_R2=np.nan,pseudo_F=F,p=p,permutations=B,lingoes_constant=lc))
  print(marker,variant,'paired n',n,'carbon partial R2',(ssr-ssf)/ssr,'p',p,flush=True)
co=pd.DataFrame(community);co['q_BH']=np.nan
for variant in co.variant.unique():
 ix=co.variant.eq(variant)&co.p.notna();co.loc[ix,'q_BH']=multipletests(co.loc[ix,'p'],method='fdr_bh')[1]
r=pd.DataFrame(reg);r['q_BH']=np.nan
for (variant,model),sub in r.groupby(['variant','model']):
 ix=sub.index[~sub.term.eq('Intercept')];r.loc[ix,'q_BH']=multipletests(r.loc[ix,'p'],method='fdr_bh')[1]
for name,d in [('community_carbon',co),('carbon_covariation',r),('cohort_summary',pd.DataFrame(summary)),('validation',pd.DataFrame(checks)),('sample_manifest',pd.concat(man))]:d.to_csv(O/(name+'.csv'),index=False)
with pd.ExcelWriter(O/'Carbon_Adjusted_Analyses.xlsx') as w:
 for f in sorted(O.glob('*.csv')):pd.read_csv(f).to_excel(w,sheet_name=f.stem[:31],index=False)
