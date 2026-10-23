from pathlib import Path
import pandas as pd,subprocess
O=Path('manuscript/draft3/carbon_paragraph_stats');a=pd.read_csv(O/'carbon_predictor_omnibus_tests.csv');b=pd.read_csv(O/'sequential_permanova_rerun.csv')
def tab(d):
 d=d.copy()
 for col in d.select_dtypes('float').columns:d[col]=d[col].map(lambda x:f'{x:.5g}')
 return '| '+' | '.join(d.columns)+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(map(str,x))+' |' for x in d.itertuples(index=False,name=None))+'\n'
paragraph='''Recorded carbon was significantly lower at 10–30 cm than at 0–10 cm. In models containing sampling date, site, and depth, carbon also differed by site (BH-adjusted q = 0.013–0.035) and sampling date (q < 0.001 for the 16S and ITS cohorts; q = 0.046 for AMF). After conflicting carbon records were excluded, depth and site remained significant in these models, whereas date did not (q = 0.198–0.417). These sensitivity results reflect both record exclusion and reduced sample/date coverage. Management and date were confounded and could not be interpreted as independent effects.

In exploratory sequential PERMANOVA models containing management, recorded carbon, and their interaction, carbon entered after management accounted for 5.2% of total variation in 16S communities (R² = 0.052, p = 0.0021, q = 0.0084), 5.8% in ITS (R² = 0.058, p = 0.0107, q = 0.0214), and 5.4% in the integrated cross-marker dataset (R² = 0.054, p = 0.0004, q = 0.0048). Carbon was not significant for AMF (R² = 0.020, p = 0.6552, q = 0.7003). The management × carbon interaction was significant for 16S (R² = 0.036, p = 0.0260, q = 0.0446), but not for ITS, AMF, or the integrated dataset (q ≥ 0.187). These models omitted sampling date, site, depth, and profile dependence and therefore do not identify an independent carbon effect. In the separate analyses conditioning on profile and depth, carbon was not significant for any individual marker (q = 0.779).'''
s='''# Astera draft 3: carbon paragraph statistics

## Corrected results paragraph

'''+paragraph+'''

## Table X. Omnibus carbon tests: date + site + depth

'''+tab(a[a.model.eq('date_site_depth')][['marker','variant','term','n','profiles','F','df_num','df_den','p','q_BH']])+'''
## Table X sensitivity. Omnibus carbon tests: management + site + depth

'''+tab(a[a.model.eq('management_site_depth')][['marker','variant','term','n','profiles','F','df_num','df_den','p','q_BH']])+'''
## Updated Table 6. Sequential community PERMANOVA

'''+tab(b[b.p.notna()][['marker','variant','term','n','R2','F','p','q_BH']])+'''
## Methods and interpretation

Carbon linear models used categorical date, site and depth, or management, site and depth. Management is aliased with date, so the two were not fitted as separately interpretable effects. Omnibus cluster-robust Wald F tests were used for each predictor, with profile-clustered standard errors, small-sample correction, and denominator degrees of freedom equal to number of profiles minus one. Date has three numerator degrees of freedom with all recorded samples and two after conflicting records are excluded. BH correction was applied across the nine predictor tests (three predictors × three marker cohorts), separately for each specification and record variant. These omnibus date q values differ from the earlier coefficient-level tests because they test date as a whole and use a different correction family. Broader field dependence remains unresolved, and profiles are assumed independent.

Community models reproduce the original formula and distances: community ~ management + carbon + management:carbon, tested sequentially using vegan::adonis2 with 9,999 unrestricted permutations, seed 20261003 reset per model. Individual-marker distances were Bray–Curtis on relative abundances. The integrated dataset comprised 42 common samples: features present in at least 5% of these samples were retained, capped at the 200 most abundant per marker. Counts received a pseudocount of one and were transformed to centered log ratios separately per marker, concatenated after explicit sample-ID alignment, and column-standardized before Euclidean distance calculation. Integrated carbon values came from the 16S metadata, as in the original analysis; they are not an independently verified cross-marker carbon measurement. Integrated feature selection and scaling were held fixed when removing conflicting samples for sensitivity analyses.

All four original sequential R² tables reproduced exactly, including the integrated dataset. The new p values differ from the historical 999-permutation estimates because they use a new fixed seed and 9,999 permutations. BH correction was added across all 12 sequential terms (three terms × four datasets), separately for each carbon record variant. Tests are correlated because marker cohorts overlap and the integrated dataset reuses them. The sequential carbon increment is conditional on earlier management in the term order, not carbon examined in isolation, nor carbon conditioned on all sampling-design factors. Full-model residual variance is the denominator of each sequential pseudo-F test. These are exploratory reproductions: unrestricted permutations do not account for profile dependence or the other omitted covariates. Prefer the separate profile- and depth-conditioned tests for that specific within-profile question; those tests do not assess the integrated dataset.

Sensitivity results changed materially: the 16S interaction no longer met q < 0.05 after conflict exclusion, ITS carbon no longer met q < 0.05, and AMF carbon became significant in the unadjusted sequential model. This reinforces the need to reconcile the chemistry records. Significance in these exploratory models does not establish causation or mediation.

The carbon type and units remain unverified, so the text uses recorded carbon instead of confirmed TOC. The quoted R² = 0.020, p = 0.658 was the historical AMF result, not an independently analyzed 18S dataset. The existing manuscript figures and files were not overwritten. Table numbering should be reconciled with the manuscript when inserting these tables.

## Reproduce

From the Astera workspace root, run source/prepare.py and source/permanova.R in this directory, then source/report.py. Python uses the saved draft3 carbon sample manifest and original count exports. R uses vegan; the software version and seed are in run_notes.txt. Source scripts, input tables, complete results, and validation are retained.
'''
(O/'Carbon_Paragraph_Statistics.md').write_text(s);(O/'Corrected_Results_Paragraph.md').write_text(paragraph)
subprocess.run(['pandoc',str(O/'Carbon_Paragraph_Statistics.md'),'-o',str(O/'Carbon_Paragraph_Statistics.docx')],check=True)
with pd.ExcelWriter(O/'Carbon_Paragraph_Statistics.xlsx') as w:
 a.to_excel(w,sheet_name='Carbon predictor tests',index=False);b.to_excel(w,sheet_name='Sequential PERMANOVA',index=False);pd.read_csv(O/'legacy_R2_validation.csv').to_excel(w,sheet_name='R2 validation',index=False)
print('Saved Word report, Excel tables, and corrected paragraph.')
