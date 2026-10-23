from common import ROOT
import pandas as pd,numpy as np
rows=[]
def compare(ref,generated,cols=None):
 a=pd.read_csv(ROOT/ref);b=pd.read_csv(ROOT/generated);assert a.shape==b.shape,(ref,a.shape,b.shape)
 assert list(a)==list(b),(ref,'columns')
 use=cols or [c for c in a if pd.api.types.is_numeric_dtype(a[c])]
 for c in use:
  aa=pd.to_numeric(a[c],errors='coerce').to_numpy(float);bb=pd.to_numeric(b[c],errors='coerce').to_numpy(float)
  passed=np.allclose(aa,bb,rtol=1e-7,atol=1e-8,equal_nan=True)
  diff=np.abs(aa-bb);delta=float(np.nanmax(diff)) if np.isfinite(diff).any() else 0.
  rows.append(dict(reference=ref,generated=generated,column=c,passed=bool(passed),max_abs_difference=delta))
for name in ['alpha_tests_rarefied','depth_selection','sample_alpha_rarefied','unique_ASV_retention_summary','per_draw_test_summary']:
 compare(f'final_reference/rarefaction/{name}.csv',f'reproduced/rarefied/rarefaction/{name}.csv')
for name in ['dbrda_model_statistics','variance_partition_summary','variance_partition_all_fractions','subset_models']:
 compare(f'final_reference/full_dbrda/{name}.csv',f'reproduced/full_dbrda/{name}.csv')
for name in ['model_statistics','permutation_diagnostics']:
 compare(f'final_reference/restricted_permanova/{name}.csv',f'reproduced/restricted_permanova/{name}.csv')
for name in ['community_carbon','carbon_covariation','cohort_summary']:
 compare(f'final_reference/adjusted_carbon/{name}.csv',f'reproduced/adjusted_carbon/{name}.csv')
a=pd.read_csv(ROOT/'final_reference/variance_without_date/variance_partition_summary.csv');b=pd.read_csv(ROOT/'reproduced/variance_without_date/variance_partition_summary.csv').set_index('marker')
col={'Management unique':'management_unique','Site unique':'site_unique','Depth unique':'depth_unique','Shared (combined)':'shared','Unexplained':'unexplained'}
for r in a.itertuples():
 value=b.loc[r.marker,col[r.component]];rows.append(dict(reference='variance_without_date',generated='variance_without_date',column=r.marker+' '+r.component,passed=bool(np.isclose(r.percent,value,atol=1e-8)),max_abs_difference=abs(r.percent-value)))
out=pd.DataFrame(rows);out.to_csv(ROOT/'validation/numerical_comparison.csv',index=False)
print('Checks:',len(out),'passed:',out.passed.sum());print(out[~out.passed].to_string(index=False));assert out.passed.all()
