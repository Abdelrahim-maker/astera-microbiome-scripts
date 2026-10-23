from pathlib import Path
import pandas as pd,numpy as np,json,re,shutil
O=Path('manuscript');R=Path('analysis_outputs');D=['16S','ITS','AMF']
read=lambda p:pd.read_csv(R/p)
def fmt(v):
 if pd.isna(v):return '—'
 if isinstance(v,(float,np.floating)):
  if v==0:return '0'
  if abs(v)<.001:return f'{v:.2e}'
  return f'{v:.4g}'
 return str(v).replace('|','/').replace('\n',' ')
def table(df):
 return '| '+' | '.join(str(c) for c in df.columns)+' |\n| '+' | '.join('---' for c in df.columns)+' |\n'+'\n'.join('| '+' | '.join(fmt(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))+'\n'
def block(label,title,df,note):return f'\n\n::: pagebreak\n:::\n\n**{label}. {title}**\n\n'+table(df)+'\n'+note+'\n'
main='''---
title: "Management-associated shifts in soil microbial composition across three amplicon datasets: carbon associations and sensitivity to unequal sampling"
author: "Authors and affiliations to be supplied"
date: "Working preprint draft | 7 September 2026"
lang: en-US
fontsize: 11pt
geometry: margin=1in
linestretch: 1.35
colorlinks: true
linkcolor: black
urlcolor: blue
bibliography: source/references.json
reference-section-title: References
---

**Draft status:** Prepared in a single-column preprint format. Field and laboratory methods require author completion. Site confounding and metadata discrepancies limit interpretation.

'''+(O/'source/manuscript_body.md').read_text()
# Main tables
sample=read('dataset_summary.csv').rename(columns={'domain':'Marker','n_samples':'Samples','n_nt':'NT','n_till':'Till','n_taxa':'Features'})
for t,label in [('NT','NT reads'),('vtill','Till reads')]:sample[label]=[pd.read_csv(O/'audit'/f'{d}_sample_alpha.csv').query('management == @t').reads.median() for d in D]
main+=block('Table 1','Analytical sample coverage and median library totals.',sample,'Feature counts refer to the selected supplied objects before downstream prevalence filters. Reads are median count-matrix totals. All selected crop labels are soybean; MF, site codes, and depth intervals require author confirmation. The integrated set contains 42 samples (33 NT, 9 till).')
a=pd.concat([read(f'{d}/{d}_alpha_wilcox.csv') for d in D]).rename(columns={'domain':'Marker','metric':'Metric','nt_median':'NT median','till_median':'Till median','p_value':'p'})
main+=block('Table 2','Alpha-diversity comparisons.',a[['Marker','Metric','NT median','Till median','p']],'Two-sided Wilcoxon rank-sum tests, normal approximation; p values are unadjusted. Sample sizes are given in Table 1. No comparison meets p < 0.05. Richness is observed ASV count; Shannon uses natural logarithms.')
rows=[];disp=read('imbalance_sensitivity/ordination_observed_permanova_permdisp.csv')
for d in D+['Joint']:
 b=read(f'{d}/{d}_nmds_permanova.csv' if d!='Joint' else 'joint_rpca_proxy_permanova.csv');v=b.loc[b.term=='Model'].iloc[0];rows.append([d,int(v.Df),int(b.loc[b.term=='Residual','Df'].iloc[0]),v['F'],v.R2,v['Pr(>F)'],disp.loc[disp.domain==d,'permdisp_p'].iloc[0] if d!='Joint' else np.nan])
main+=block('Table 3','Primary community-level management tests.',pd.DataFrame(rows,columns=['Marker','df','Residual df','Pseudo-F','R²','p','Dispersion p']),'PERMANOVA uses 999 unrestricted permutations. Marker-specific tests use Bray–Curtis dissimilarities of relative abundance; Joint uses Euclidean distance on standardized concatenated CLR features. Dispersion p values come from the sensitivity run. Primary CSVs are the authority for the reported main-model p values; archived reruns may differ slightly because permutation draws differ.')
rows=[]
for d in D:
 a=read(f'{d}/{d}_maaslin2_nt_vs_till_results.csv');s=a[a.qval<=.25];rows.append([d,len(a),len(s),(a.qval<=.05).sum(),(s.coef<0).sum(),(s.coef>0).sum()])
main+=block('Table 4','Treatment-only MaAsLin2 discovery counts.',pd.DataFrame(rows,columns=['Marker','Tested','q ≤ 0.25','q ≤ 0.05','Higher NT¹','Higher till¹']),'¹ Direction counts use q ≤ 0.25. BH correction was applied to the feature tests in each assay. q ≤ 0.25 is exploratory. Features and relative-abundance coefficient estimates are provided in Table S3 and the complete data workbook.')
a=read('network_stats_summary.csv').rename(columns={'label':'Network','nodes':'Nodes','edges':'Edges','density':'Density','mean_degree':'Mean degree','transitivity':'Transitivity','modularity':'Modularity'})
main+=block('Table 5','Observed association-network summaries.',a,'Networks retain Spearman |ρ| ≥ 0.6 and BH q < 0.05 after count-plus-one CLR transformation and feature filtering. Nodes exclude isolates. Modularity uses absolute correlation weights. These summaries do not measure direct interactions, stability, or field replication.')
a=read('carbon_analysis/carbon_permanova_by_domain.csv');j=read('carbon_analysis/carbon_joint_permanova.csv');j['domain']='Joint';a=pd.concat([a,j]);a=a[(a.test_type=='sequential') & a.term.isin(['trt2','C','trt2:C'])].copy();a.term=a.term.map({'trt2':'Management','C':'C','trt2:C':'Management × C'});a=a[['domain','term','Df','F','R2','Pr(>F)']];a.columns=['Marker','Term','df','Pseudo-F','R²','p']
main+=block('Table 6','Sequential carbon-inclusive PERMANOVA.',a,'Terms were entered in the displayed order; the management term is not adjusted for C. Residual degrees of freedom are 40, 39, 40, and 38 for 16S, ITS, AMF, and Joint. C units are unspecified. Models omit site and depth. AMF chemistry differs for 20 matched samples and requires reconciliation.')
# Main figures: use original bitmaps unchanged where applicable.
figs=json.loads((O/'source/figure_index.json').read_text());lookup={f['source']:f for f in figs}
def orig(p):return lookup['analysis_outputs/'+p]['file']
F=[
('figures/figure1_alpha.png','Alpha diversity in the supplied count matrices.','Observed richness (A–C) and Shannon entropy (D–F), ordered 16S, ITS, and AMF. Each point is an analytical sample; NT and till counts are 33/11, 33/10, and 34/10. Boxes show medians and interquartile ranges; whiskers extend to 1.5 times the interquartile range. Jitter is cosmetic. p values are the original two-sided Wilcoxon results, not corrected for multiple comparisons. Values were recalculated from the supplied counts and matched the source medians.'),
(orig('16S/16S_nmds.png'),'16S composition under NT and till.','Original Bray–Curtis NMDS from relative abundances; points are 33 NT and 11 till observations. Shaded ellipses are the plotting function’s default 95% t-type group ellipses and are not a test of treatment separation. Stress = 0.141. The current source CSV gives PERMANOVA R² = 0.04695, p = 0.006 with 999 unrestricted permutations. Site and depth were not included. ITS and AMF ordinations are reproduced in the supplementary atlas.'),
(orig('joint_rpca_proxy_plot.png'),'Integrated ordination of shared samples.','Original PCA of standardized, concatenated count-plus-one CLR marker features for 42 shared samples (33 NT, 9 till). Despite the legacy “RPCA proxy” label, this is conventional standardized CLR-PCA, not robust Aitchison PCA. Axis percentages are PCA variance proportions. PERMANOVA on the standardized matrix gave R² = 0.05095, p = 0.004. The two fungal assays are not independent taxonomic domains.'),
('figures/figure4_differential.png','Feature discovery depends on the adjusted-p threshold.','Treatment-only MaAsLin2 counts for 16S, ITS, and AMF. Dark bars show q ≤ 0.05; pale segments show 0.05 < q ≤ 0.25. Total exploratory counts are 34, 43, and 3 from 915, 127, and 87 tested features, respectively. Bars show counts, not effect sizes, and therefore have no sampling-error bars. Full coefficients, standard errors, p values, and q values are in the supplement.'),
('figures/figure5_networks.png','Observed network topology by management.','Connected nodes, mean degree, and density for the three marker datasets and the cross-marker (Joint) graph. Green denotes NT and orange till. Values are descriptive summaries of one inferred graph per group; error bars are not available. Candidate feature sets differ between groups, isolated features are excluded, and sample imbalance affects edge detection. The large observed AMF contrast contracts in balanced subsamples (Table S7).'),
('figures/figure6_sensitivity.png','Community-level sensitivity to balanced subsampling.','Each point represents one of 20 balanced subsets; boxes show medians and interquartile ranges and whiskers extend to 1.5 times that range. Left: PERMANOVA R². Right: permutation p values, with a dashed line at 0.05. Group sizes are 11/11 for 16S, 10/10 for ITS and AMF, and 9/9 for Joint. The fractions below 0.05 are 18/20, 7/20, 4/20, and 20/20. Subsets overlap and are not independent biological replications. Site confounding remains.'),
(orig('carbon_analysis/carbon_permanova_r2_heatmap.png'),'Carbon-inclusive community models.','Original visualization of sequential PERMANOVA terms from the supplied carbon analysis. Values and significance annotations refer to management, C, and management × C in that order; Table 6 provides numeric results. Nominal term-level p values were not corrected across the family of community models. C units are unspecified, site and depth are omitted, and AMF chemistry is unresolved.'),
(orig('ITS/ITS_network_till_carbon_overlay.png'),'ITS till network with carbon-associated feature overlays.','Original exploratory graph for the ten ITS till observations. Borders identify C associations (green), interaction associations (orange), or both (purple) in carbon-inclusive models at q ≤ 0.25. Fill colors and square shapes indicate treatment-only associations as defined in the source legend; larger nodes identify degree-based hubs. No ITS interaction term met q ≤ 0.25. Across both ITS networks, eight of nine C-associated ASVs occur as nodes, one is a hub, and three are treatment-associated. Border associations come from the full carbon model, not a till-only carbon-slope test. Neither edges nor hubs establish ecological interactions.')]
for i,(path,title,cap) in enumerate(F,1):
 main+=f'\n\n::: pagebreak\n:::\n\n**Figure {i}. {title}**\n\n![]({path}){{width=6.3in height=5.3in}}\n\n{cap}\n'
main+='\n\n::: pagebreak\n:::\n\n## References\n\n::: {#refs}\n:::\n'
# Figure atlas references are explicit and exhaustive.
indexpara='\n\nThe supplementary figure atlas reproduces every original image as Figures S1–S'+str(len(figs))+'. The figure index maps each number to its exact source file. The complete statistical output is provided as data tables D001–D'+f'{len(json.loads((O/"source/data_index.json").read_text())):03d}'+', indexed in the workbook and supplement. Figures 1, 4–6 are new visual summaries of supplied data; all other displayed source figures are unchanged.\n'
main=main.replace('## Data and code availability',indexpara+'\n## Data and code availability')
(O/'manuscript.md').write_text(main)
# Supplement with numbered summary tables plus complete atlas.
supp='''---
title: "Supplementary material: Management-associated shifts in soil microbial composition"
author: "Working preprint draft"
date: "7 September 2026"
fontsize: 10pt
geometry: margin=0.8in
linestretch: 1.1
---

This supplement preserves the supplied analyses and records manuscript quality-control findings. Tables S1–S12 summarize design, features, sensitivity, and carbon models. Data tables D001 onward contain every original CSV export, including full model uncertainty fields and all bootstrap records. The original figures are reproduced unchanged; legacy titles, clipping, and labels are retained for provenance, with interpretive qualifications in captions. The manuscript presents associations rather than causal treatment effects.
'''
rows=[]
for d in D:
 m=pd.read_csv(O/'audit'/f'{d}_metadata.csv');m=m[(m.treatment=='MF')&m.management.isin(['NT','vtill'])]
 for field in ['site','depth']:
  for level,g in m.groupby(field):rows.append([d,field,level,(g.management=='NT').sum(),(g.management=='vtill').sum()])
supp+=block('Table S1','Sample distribution across site and depth codes.',pd.DataFrame(rows,columns=['Marker','Field','Code','NT','Till']),'All selected crop labels are soybean. Site and depth codes are not expanded into unverified physical definitions. Management and site are partially confounded.')
a=pd.read_csv(O/'audit/metadata_discrepancies.csv');summ=a.groupby(['comparison','field']).size().reset_index(name='Disagreements');supp+=block('Table S2','Cross-marker metadata discrepancies.',summ,'Matched samples were compared by exact row identifier; numeric fields used floating-point tolerance. Values differ for 20 AMF C, N, and pH records. Full sample-level discrepancies are in audit/metadata_discrepancies.csv. The correct source values have not been determined.')
rows=[]
for d in D:
 a=read(f'{d}/{d}_maaslin2_nt_vs_till_results.csv').sort_values('qval').head(5)
 for _,r in a.iterrows():rows.append([d,r.feature,r.coef,r.stderr,r.pval,r.qval])
supp+=block('Table S3','Five lowest-q treatment features per marker.',pd.DataFrame(rows,columns=['Marker','ASV','Coefficient','SE','p','q']),'Selection is by q value, irrespective of whether the fifth feature meets a threshold. Positive coefficients indicate higher transformed relative abundance in till. This selected table does not replace all 1,129 feature rows and their taxonomy in the complete workbook.')
a=read('enriched_hub_fisher_tests.csv');supp+=block('Table S4','Archived treatment-feature and hub overlap tests.',a,'These source-table universes are restricted to the network-related overlap analysis and therefore differ from the full MaAsLin2 feature counts. The p values are archived one-sided Fisher enrichment outputs. No overrepresentation was detected; source code for generation of some overlap-only exports was not supplied.')
a=read('imbalance_sensitivity/ordination_balanced_summary.csv')[['domain','median_permanova_r2','median_permanova_p','frac_permanova_p_lt_0_05']];a.columns=['Marker','Median R²','Median p','Fraction p < 0.05'];j=read('imbalance_sensitivity/joint_balanced_summary.csv').iloc[0];a.loc[len(a)]=['Joint',j.median_r2,j.median_p,j.frac_p_lt_0_05];supp+=block('Table S5','Balanced community-level sensitivity.',a,'Twenty overlapping balanced subsets per analysis. Fractions describe the supplied runs; they are not confidence levels or power estimates.')
a=read('imbalance_sensitivity/maaslin_balanced_summary.csv');a.columns=['Marker','Median discoveries','Minimum','Maximum'];supp+=block('Table S6','Balanced MaAsLin2 discovery counts.',a,'Discovery threshold q ≤ 0.25. Full iteration records and feature recurrence are included as indexed data tables. A fractional median is possible because there are 20 iterations.')
a=read('imbalance_sensitivity/network_balanced_summary.csv')[['domain','label','median_nodes','median_edges','median_mean_degree','median_density']];a.columns=['Marker','Group','Median nodes','Median edges','Median degree','Median density'];supp+=block('Table S7','Balanced network topology.',a,'The smaller till group is retained in every iteration. Node-degree p values are preserved in the complete tables but are not valid independent field-level comparisons. Sample-size changes also alter prevalence and correlation-threshold behavior.')
a=read('carbon_analysis/carbon_treatment_wilcox.csv')[['domain','nt_median','till_median','p_value']];a.columns=['Marker','NT median C','Till median C','p'];supp+=block('Table S8','Carbon group comparisons.',a,'C measurement and units are unverified. These assay sample sets overlap. AMF chemistry differs from the other assay metadata for many matched samples; values are reproduced without choosing a replacement.')
a=read('carbon_analysis/carbon_maaslin_effect_summary.csv');a.effect_type=a.effect_type.replace({'carbon_main_NT_slope':'C slope (NT)','carbon_slope_difference_till_vs_NT':'Slope interaction','treatment_main':'Treatment at C = 0'});a.columns=['Marker','Term','Tested','q ≤ 0.25'];supp+=block('Table S9','Carbon-inclusive feature discovery.',a,'Continuous predictors C and C_by_till were standardized separately. Interaction coefficients require rescaling before combination with the NT slope. q values belong to the carbon-inclusive output and should not be confused with treatment-only model q values.')
a=read('carbon_analysis/carbon_taxa_network_enrichment_overlap_summary.csv');a=a[['domain','effect_type','n_carbon_sig','n_in_network','n_in_hubs','n_enriched_in_network']];a.effect_type=a.effect_type.replace({'carbon_main_NT_slope':'C slope','carbon_slope_difference_till_vs_NT':'Interaction'});a.columns=['Marker','Term','Associated','In graph','Hubs','Treatment overlap'];supp+=block('Table S10','Carbon-associated features and graph overlap.',a,'Counts are term-specific; C and interaction categories can overlap within an assay. AMF has no q ≤ 0.25 carbon terms. ITS graph membership is the union across treatment networks.')
a=read('carbon_analysis/carbon_variance_partition_by_domain.csv');a.columns=['Marker','Fraction','Adjusted R²'];supp+=block('Table S11','Additive Hellinger variance partitioning.',a,'[a] is management unique of C, [b] is C unique of management, [c] is shared, and [d] is residual. Fractions are descriptive, can be negative after adjustment, and are not interaction PERMANOVA effects.')
a=read('carbon_analysis/network_metric_carbon_associations.csv')[['domain','trt2','cor_mean_degree','q_mean_degree','cor_modularity','q_modularity']];a.columns=['Marker','Group','Degree r','Degree q','Modularity r','Modularity q'];supp+=block('Table S12','Archived bootstrap network–carbon correlations.',a,'Each group has 80 resamples of up to ten observations with replacement. Bootstrap samples overlap and are not independent field replicates. These adjusted values should not be interpreted as confirmatory treatment-level inference. Full density results and unadjusted p values are preserved in the data workbook.')
supp+='\n\n::: pagebreak\n:::\n\n## Complete statistical-table index\n\nEach D-number corresponds to one worksheet in supplementary_tables.xlsx. All rows and columns are retained; CSV/TSV originals are also in supplementary_data/. File paths below are relative to analysis_outputs/.\n\n'
for entry in json.loads((O/'source/data_index.json').read_text()):supp+=f"- **{entry['id']}**: `{entry['source'].replace('analysis_outputs/','')}` ({entry['rows']} rows, {entry['columns']} columns).\n"
def caption(f):
 p=f['source'];name=Path(p).stem;d=next((x for x in D if '/'+x+'/' in p),'Combined');base=f'Original source: `{p}`. '
 if 'network' in name and 'bootstrap' not in name and 'overlap' not in name and 'dotplot' not in name:
  base+='Exploratory association network or network panel. Edges require Spearman |ρ| ≥ 0.6 and BH q < 0.05 after count-plus-one CLR transformation. Nodes exclude isolated candidate features. Hubs follow the source degree rule and are not validated keystone taxa. '
  if 'carbon_overlay' in name:base+='Borders encode carbon-model q ≤ 0.25 associations (green: NT C slope; orange: interaction; purple: both). Fills/shapes encode treatment-only enrichment. Associations are from the full model, not a treatment-specific correlation test. '
  elif 'overlay' in name:base+='Directional outlines encode treatment-only MaAsLin2 q ≤ 0.25 (green: higher NT; orange: higher till). Common/important-node markers in panel versions follow the embedded legend. '
  else:base+='Blue/red edges identify positive/negative correlations where shown. '
  base+='Topology is sensitive to sample count and filtering; site and depth are unadjusted.'
 elif 'nmds' in name:base+='Bray–Curtis NMDS on relative abundance. Points are analytical samples; shaded ellipses are default 95% t-type group ellipses. Titles report stress and 999-permutation management-only PERMANOVA. The combined legacy panel has cropped title margins; individual panels are supplied separately. Nonsignificant dispersion tests do not remove site confounding.'
 elif 'rarefaction' in name:base+='Descriptive rarefaction curves from the processed count tables; NT is green and till orange. These do not demonstrate that all soil diversity was sampled. The combined image is a legacy panel of the individual curves.'
 elif 'shannon' in name or 'observed' in name or 'alpha_diversity' in name:base+='Original alpha-diversity display. Boxes indicate medians/interquartile ranges with default 1.5-IQR whiskers; points represent samples. Wilcoxon results are in Table 2. Source objects produced singleton warnings. Assay counts are 33/11 (16S), 33/10 (ITS), and 34/10 (AMF), NT/till.'
 elif 'lollipop' in name:base+='Treatment-only MaAsLin2 coefficients for up to 30 q ≤ 0.25 features. Positive values indicate higher transformed relative abundance in till; negative values indicate higher NT. Taxonomy labels are supplied annotations. These plots do not display all tested features or coefficient uncertainty; full standard errors and adjusted values are in the workbook.'
 elif 'joint_rpca' in name:base+='PCA of standardized concatenated CLR features for 42 shared observations (33 NT, 9 till). The legacy RPCA label denotes a proxy only; robust Aitchison matrix completion was not used.'
 elif 'carbon_by_treatment' in name:base+='Original C distributions by group. Measurement method and units are unknown. Wilcoxon comparisons appear in Table S8; overlapping assay samples are not independent cohorts.'
 elif 'permanova' in name:base+='Sequential carbon PERMANOVA summary; management, C, and interaction were entered in that order. Nominal term p values are exploratory; Table 6 gives the exact statistics.'
 elif 'effect_counts' in name:base+='Carbon-inclusive MaAsLin2 discovery counts at q ≤ 0.25. C denotes the NT slope and C_by_till a scaled difference in slopes. Separately standardized predictors prevent naive coefficient addition.'
 elif 'bootstrap' in name:base+='Network mean degree plotted against mean C for 80 overlapping bootstrap resamples per marker–treatment group, each of up to ten draws with replacement. Fitted trends and any bands are descriptive; resamples are not independent field observations.'
 elif 'correlations' in name:base+='Selected feature–carbon Pearson correlations among hub/enriched features. Selection and reuse of the same samples make these exploratory; full coefficients, p values, and adjusted values are archived.'
 elif 'overlap' in name:base+='Carbon q ≤ 0.25 feature overlap with graph membership, hubs, and treatment association. Counts are term-specific, not necessarily distinct ASVs across terms. See Table S10.'
 elif 'dotplot' in name:base+='Ranked ITS carbon-associated network features. Supplied annotations and graph degrees identify follow-up candidates, not verified ecological functions. Eight ITS features meet the overlap definition.'
 else:base+='Original analytical visualization; supporting numeric output is in the complete workbook.'
 if 'carbon' in p:base+=' Carbon outputs omit site/depth adjustment; AMF chemistry has unresolved matched-sample discrepancies.'
 return base
supp+='\n\n## Complete figure atlas\n'
for f in figs:
 title=Path(f['source']).stem.replace('_',' ')
 supp+=f"\n\n::: pagebreak\n:::\n\n**Figure S{f['number']}. {title}.**\n\n![]({f['file']}){{width=6.6in height=5.6in}}\n\n{caption(f)}\n"
(O/'supplementary_material.md').write_text(supp)
# A compact mapping is convenient outside the long atlas.
pd.DataFrame(figs).to_csv(O/'figure_index.csv',index=False)
for p in Path('.').glob('*.R'):shutil.copy2(p,O/'source'/p.name)
print('Manuscript and supplement assembled.')
