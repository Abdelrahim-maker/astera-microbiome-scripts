from pathlib import Path
import json
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
ROOT=Path(__file__).resolve().parents[1];T=ROOT/'manuscript/tables';checks=[]
def check(name,ok):
 assert ok,name
 checks.append(dict(check=name,passed=True))
for name in ['alpha_tests_rarefied.csv','unique_ASV_retention_summary.csv']:
 a=pd.read_csv(ROOT/'reproduced/rarefied/rarefaction'/name)
 b=pd.read_csv(ROOT/'inputs/final_reference/rarefaction'/name)
 check(name+' columns',list(a.columns)==list(b.columns))
 for c in a:
  check(name+' '+c,np.allclose(a[c],b[c],rtol=1e-9,atol=1e-10,equal_nan=True) if pd.api.types.is_numeric_dtype(a[c]) else a[c].equals(b[c]))
x=pd.read_csv(T/'network_balanced_subsamples_stats.csv');s=pd.read_csv(T/'S6_network_balanced.csv')
for row in s.itertuples():
 z=x[x.domain.eq(row.domain)&x.label.eq(row.label)]
 check(row.domain+' '+row.label+' 20 draws',len(z)==20 and z['iter'].nunique()==20)
 for c in ['nodes','edges','density','mean_degree']:check(row.domain+' '+row.label+' median '+c,np.isclose(z[c].median(),getattr(row,'median_'+c)))
b=pd.read_csv(T/'ordination_balanced_subsamples_permanova_permdisp.csv');s=pd.read_csv(T/'S7_beta_balanced.csv')
for row in s.itertuples():
 z=b[b.domain.eq(row.domain)]
 check(row.domain+' beta 20 draws',len(z)==20)
 check(row.domain+' beta frequency',np.isclose((z.permanova_p<.05).mean(),row.frac_permanova_p_lt_0_05))
 check(row.domain+' beta median R2',np.isclose(z.permanova_r2.median(),row.median_permanova_r2))
x=pd.read_csv(T/'MaAsLin2_all_results.csv')
for marker in ['16S','ITS','AMF']:
 z=x[x.marker.eq(marker)]
 check(marker+' MaAsLin2 BH',np.allclose(z.qval,multipletests(z.pval,method='fdr_bh')[1]))
(ROOT/'validation').mkdir(exist_ok=True)
(ROOT/'validation/numerical_validation.json').write_text(json.dumps({'checks':checks,'scope':'Rarefaction reproduction and retained-table checks. PERMANOVA and network reconstruction assert agreement in their own scripts. R MaAsLin2 and historical balanced subsets were not refitted.'},indent=2))
print('PASS:',len(checks),'retained-result checks')
