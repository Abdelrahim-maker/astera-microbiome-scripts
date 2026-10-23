from pathlib import Path
import pandas as pd
import subprocess
O=Path('manuscript/draft3/carbon_adjusted')
c=pd.read_csv(O/'community_carbon.csv');r=pd.read_csv(O/'carbon_covariation.csv')
def table(d):
 return '| '+' | '.join(d.columns)+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in d.itertuples(index=False,name=None))+'\n'
a=c[c.model.eq('within_profile_carbon_adjusted_depth')].copy();a['partial_R2']=a.partial_R2.round(4);a=a[['marker','variant','n','profiles','partial_R2','p','q_BH']]
b=r[r.model.eq('management_depth_site') & r.term.str.contains('management|Depth')].copy();b['term']=b.term.map(lambda x:'Till minus NT' if 'management' in x else '10–30 minus 0–10 cm');b=b[['marker','variant','term','estimate','p','q_BH']];b['estimate']=b.estimate.round(4);b['p']=b.p.map(lambda x:f'{x:.4g}');b['q_BH']=b.q_BH.map(lambda x:f'{x:.4g}')
paragraph='''The recorded carbon variable was associated with soil depth, whereas its association with management was sensitive to the carbon records used and was not significant after conflicting records were excluded. In partial community analyses controlling for soil profile and depth, carbon was not significantly associated with community composition in 16S, ITS, or AMF (BH-adjusted q = 0.779 for each marker). Thus, these analyses do not establish carbon as an explanation for the management-associated community patterns. Management was confounded with sampling date, and the carbon measurement type and units remain unverified. More generally, tillage can jointly alter carbon distribution and microbial communities across soil depths, complicating attribution to individual factors [1,2].'''
text='''# Astera draft 3: adjusted carbon analyses

## Answer and scope

The earlier analyses tested management, recorded carbon, and their interaction, but omitted depth, site, sampling date, and profile dependence. Unadjusted management comparisons were nonsignificant for all three markers. The new analyses below explicitly evaluate carbon covariation and within-profile community associations. They do not demonstrate mediation or causal explanation of management effects.

## Replacement discussion paragraph

'''+paragraph+'''

## Results

Recorded carbon was lower at 10–30 cm than at 0–10 cm in all marker-specific cohorts after accounting for management and the original two-level site code, with profile-clustered standard errors. The depth association persisted when date replaced management and after excluding conflicting carbon records. In the management + depth + site models, management met BH q < 0.05 only with the recorded AMF carbon values (q = 0.0268); the corresponding q values were 0.0636 for 16S and 0.0885 for ITS. After exclusion of conflicting records, none met q < 0.05 (16S 0.215; ITS 0.285; AMF 0.0855). These partly overlapping marker cohorts are sensitivity comparisons, not three independent carbon datasets.

The within-profile community tests used 20 complete profiles for 16S, 19 for ITS, and 20 for AMF. Carbon accounted for 7.1%, 7.7%, and 4.6% of the community variation remaining after profile and depth, respectively; none was significant. These partial R² values are fractions of residual variation, not total community variation. Exclusion of conflicting records also yielded no significant associations. Failure to detect an association does not establish absence of a carbon effect, particularly with the smaller sensitivity cohorts.

## Table 1. Within-profile carbon association with community composition

'''+table(a)+'''
## Table 2. Carbon covariation with management and depth

Estimates are in the recorded, unverified carbon units. BH correction includes all non-intercept coefficients across the three markers within each model specification and record variant.

'''+table(b)+'''
## Methods and reproducibility

The original draft 3 cohorts were retained: 44 16S, 43 ITS, and 44 AMF samples. Carbon was joined by soil UID using the existing sample-table values separately for each marker; no carbon source was declared authoritative. Workbook Site_Name defined profiles, and original site codes defined the two-level site covariate. The sensitivity analysis excluded every soil UID flagged as having conflicting carbon records, reducing the full cohorts to 24, 23, and 24 samples. The paired sensitivity subsets contained 11, 10, and 11 complete profiles.

Carbon covariation was assessed with ordinary least squares: recorded C ~ management + categorical depth + site, and recorded C ~ categorical date + categorical depth + site. Standard errors were clustered by profile, using the small-sample correction and t inference with number of profiles minus one degrees of freedom. These exploratory tests assume profiles are independent clusters; broader field dependence remains unresolved. Date and management cannot be interpreted as independently estimated predictors because tillage is aliased with sampling date.

Community counts were aligned by sample ID, converted to relative abundance, and compared with Bray–Curtis distances. Lingoes correction was applied when needed to obtain Euclidean principal coordinates. Descriptive full-cohort models included management, date, site, depth, and carbon; each unique increment was computed after all other terms. Management added zero degrees of freedom after date. No p values were assigned to these descriptive full-cohort increments.

The inferential community model was community ~ carbon + Condition(profile + categorical depth), restricted to complete two-depth profiles. Profile absorbs all between-profile differences, including management, date, and site, which are constant within each profile. The test therefore concerns within-profile carbon variation after the common depth effect; it cannot identify a between-profile carbon effect or carbon mediation of management. Reduced-model residuals were swapped between depths within profiles and added back to fitted values in 9,999 Monte Carlo permutations (seed 20261003). Models were refitted for each permutation, and p = (exceedances + 1)/10,000. BH correction was across three markers separately for each record variant. These residual permutation tests are approximate conditional tests requiring residual exchangeability within profiles; profile groups must also be independent.

Python checks verified corrected distances, projection matrices, least-squares residual sums, and partial R². All six paired model F statistics and partial R² values were independently checked with R/vegan. Scripts, sample manifest, full numeric tables, and validation outputs are saved alongside this report. Reproduce from the Astera root with python manuscript/draft3/carbon_adjusted/source/analyze.py, then Rscript manuscript/draft3/carbon_adjusted/source/validate.R and the report builder. Original draft files were preserved.

## References

1. Effects of transitioning from conventional to organic farming on soil organic carbon and microbial community: a comparison of long-term non-inversion minimum tillage and conventional tillage. Biology and Fertility of Soils. https://doi.org/10.1007/s00374-024-01796-y

2. Conservation tillage facilitates the accumulation of soil organic carbon fractions by affecting the microbial community in an eolian sandy soil. Frontiers in Microbiology. https://doi.org/10.3389/fmicb.2024.1394179

3. vegan: permutation tests for constrained ordination. https://vegandevs.github.io/vegan/reference/anova.cca.html

References 1–2 support the general observation that tillage, depth, carbon, and microbial communities can change together; they do not validate the specific carbon records or establish confounding in this dataset.
'''
assert pd.read_csv(O/'vegan_validation.csv').passed.all()
(O/'Carbon_Adjusted_Report.md').write_text(text)
(O/'Replacement_Paragraph.md').write_text(paragraph)
subprocess.run(['pandoc',str(O/'Carbon_Adjusted_Report.md'),'-o',str(O/'Carbon_Adjusted_Report.docx')],check=True)
with pd.ExcelWriter(O/'Carbon_Adjusted_Analyses.xlsx') as w:
 for name in ['community_carbon','carbon_covariation','cohort_summary','validation','vegan_validation','sample_manifest']:pd.read_csv(O/(name+'.csv')).to_excel(w,sheet_name=name,index=False)
print('Saved Word report, replacement paragraph, Excel tables, and reproducible analyses.')
