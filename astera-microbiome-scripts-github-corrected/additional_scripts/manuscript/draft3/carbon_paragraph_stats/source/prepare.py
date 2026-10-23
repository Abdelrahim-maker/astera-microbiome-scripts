from pathlib import Path
import pandas as pd,numpy as np,json
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
O=Path('manuscript/draft3/carbon_paragraph_stats');m=pd.read_csv('manuscript/draft3/carbon_adjusted/sample_manifest.csv');rows=[]
for marker in ['16S','ITS','AMF']:
 d=m[m.marker.eq(marker)].copy();info=json.loads(Path(f'manuscript/audit/{marker}_matrix_metadata.json').read_text());a=np.load(f'manuscript/audit/{marker}_counts.npz')['values'];ids=[(i,s) for i,s in enumerate(info['dimnames'][0]) if s in set(d.sample_id)];d=d.set_index('sample_id').loc[[s for i,s in ids]].reset_index();pd.DataFrame(a[[i for i,s in ids]],index=d.sample_id).to_csv(O/f'{marker}_counts.csv');d.to_csv(O/f'{marker}_metadata.csv',index=False)
 for variant in ['recorded','exclude_conflicting']:
  x=d if variant=='recorded' else d[~d.carbon_conflict]
  for label,formula in [('date_site_depth','carbon ~ C(sample_date) + C(site) + C(Depth_cm)'),('management_site_depth','carbon ~ C(management_group) + C(site) + C(Depth_cm)')]:
   f=smf.ols(formula,x).fit(cov_type='cluster',cov_kwds={'groups':x.Site_Name,'use_correction':True},use_t=True)
   for term in ['sample_date' if label.startswith('date') else 'management_group','site','Depth_cm']:
    ix=[i for i,t in enumerate(f.params.index) if t.startswith('C('+term+')')];L=np.eye(len(f.params))[ix];w=f.f_test(L)
    rows.append(dict(marker=marker,variant=variant,model=label,term=term,n=len(x),profiles=x.Site_Name.nunique(),F=float(w.fvalue),df_num=float(w.df_num),df_den=float(w.df_denom),p=float(w.pvalue),model_R2=f.rsquared))
r=pd.DataFrame(rows);r['q_BH']=np.nan
for _,g in r.groupby(['variant','model']):r.loc[g.index,'q_BH']=multipletests(g.p,method='fdr_bh')[1]
r.to_csv(O/'carbon_predictor_omnibus_tests.csv',index=False)
