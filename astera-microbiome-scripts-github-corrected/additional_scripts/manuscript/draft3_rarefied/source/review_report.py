from pathlib import Path
import zipfile,lxml.etree as E
import pandas as pd,numpy as np
from docx import Document
from docx.shared import Inches,Pt
O=Path('manuscript/draft3_rarefied'); ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
responses={}
def add(ids,status,response):
 for i in ids:responses[i]=(status,response)
add(range(5),'Administrative','Removed draft email and scheduling discussion from manuscript; produced Word draft. No coauthor messages were sent.')
add([5,41,42,43,44,45,46],'Addressed','Added date-inclusive dbRDA and four-set adjusted variance partitioning (main Figures 3–4; Tables 7–8), date-colored PCoA (Figure S54), and date diagnostics (Table S13). Added Michigan temporal-variation literature. Corrected season claim: these samples span ten days in late autumn, not multiple seasons. Management cannot be separated from 28 November date.')
add([6,19],'Addressed','Rebuilt native editable Word tables with shorter headings, adjusted column widths, appropriate precision, and repeated header rows.')
add([7,9,11,13],'Partly addressed; source data missing','Pooled read medians across NT/till per marker and supplied sample-level retained count totals plus date/depth metadata. Pre-QC and intermediate-stage sequencing read counts are unavailable; explicitly marked missing, never inferred from retained matrix totals.')
add([8,10,12,14],'Addressed','Read counts use comma-separated decimal numbers rather than scientific notation. The pooled 18S median ends in .5 because the sample count is even.')
add([15,16,17],'Addressed','Defined retained reads and ASVs; moved cohort information into Methods. Distinguished observed ASVs from zero-only matrix columns. Clarified that retained totals are not raw/pre-QC yields.')
add([18],'Addressed','Added BH q values across all six alpha-diversity comparisons, retaining raw p values. Replaced the prior unrarefied alpha comparison with 200-draw rarefied alpha diversity at explicitly selected depths; ASVs are not clustered OTUs.')
add([20],'Addressed with precise terminology','Table 4 now describes differentially associated ASVs and points to full coefficient outputs. ASVs are not necessarily distinct taxonomic lineages; exploratory q ≤ 0.25 is explicit.')
add(range(21,30),'Addressed','Network counts are integers; density, transitivity, and modularity use three decimals, mean degree two. Full precision retained in source tables.')
add([30,57],'Addressed','Moved graph construction and inferential limitations to Methods/Discussion; simplified figure and table captions while retaining definitions needed to read them.')
add([31,32],'Addressed','Rebuilt main figures and tables with 16S, ITS, and 18S (AMF) labels. Defined AMF as the biological target of the 18S assay; biological AMF discussion retains that term.')
add([33,58],'Addressed','Rebuilt boxplots with black median lines; treatment colors are consistent.')
add([34],'Addressed','Rewrote alpha-diversity caption to explain points, whiskers, cosmetic horizontal offsets, raw p/BH q, sample-level rarefied means, and varying scales.')
add([35,36,37,47],'Addressed','Figure 2 now includes 16S, ITS, and 18S PCoA panels plus cross-marker standardized CLR-PCA as panel D.')
add([38,48,50,52],'Addressed','Removed MF and obsolete RPCA proxy wording from revised main displays; used fresh-frozen in cohort description. Removed unnecessary Original language from revised captions. Archival supplement is explicitly distinguished from revised main figures.')
add([39],'Addressed without mislabeling uncertainty','Removed group ellipses from the revised ordination. Sample concentration ellipses are not automatically confidence intervals for group centroids.')
add([40],'Addressed with design qualification','Added site/depth/date-inclusive dbRDA and explicit permutation assumptions. Management-only PERMANOVA remains labeled unadjusted; adding date makes unique management non-identifiable.')
add([49,51],'Addressed','Specified each transformation by analysis. Explained why PCA on log-ratios is a descriptive option; feature standardization changes ordinary Aitchison geometry. Cross-marker R² independently reconstructed from count exports. Main marker displays use Bray–Curtis PCoA.')
add([53,54],'Addressed','Axis variance percentages placed in axis labels; model statistics placed above panels. Captions distinguish corrected PCoA geometry from the uncorrected distances used by saved PERMANOVA tests.')
add([55],'Addressed','Replaced stacked counts with narrow, side-by-side bars labeled as totals for q ≤ 0.05 and q ≤ 0.25; explicitly explained inclusion.')
add([56],'Addressed','Standardized new narrative, figures, and tables to cross-marker. Historical filenames remain unchanged for provenance.')
add([59],'Figure clarified; discussion remains for authors','Redesigned sensitivity display and explained overlapping subsets, fixed till observations, detection counts, and why frequencies are not power or independent replications. The requested lab discussion cannot be performed through document editing.')
add([60],'Addressed','Rebuilt carbon plot and caption as sequential PERMANOVA. Removed unsupported total organic carbon label and any implication of carbon-adjusted management effects. Missing carbon units/fraction and discordant AMF chemistry remain explicit.')
assert set(responses)==set(range(61))
with zipfile.ZipFile('analysis_outputs/manuscript_dbrda/Astera_manuscript_draft2_comments(1).docx') as z:cs=E.fromstring(z.read('word/comments.xml'))
rows=[]
for c in cs:
 i=int(c.get('{'+ns['w']+'}id'));rows.append({'comment_id':i,'comment':''.join(c.xpath('.//w:t/text()',namespaces=ns)),'status':responses[i][0],'response':responses[i][1]})
a=pd.DataFrame(rows).sort_values('comment_id');a.to_csv(O/'audit/comment_response_log.csv',index=False)
checks=[]
for fn in ['dbrda_model_statistics','variance_partition_summary','variance_partition_all_fractions','subset_models']:
 old=pd.read_csv(f'analysis_outputs/manuscript_dbrda/{fn}.csv');new=pd.read_csv(O/f'audit/reproduced_dbrda/{fn}.csv');assert list(old.columns)==list(new.columns) and old.shape==new.shape
 for col in old.select_dtypes('number'):assert np.allclose(old[col],new[col],equal_nan=True,atol=1e-12)
 checks.append(f'{fn}: all numeric values reproduced within 1e-12.')
report=Document();report.add_heading('Draft 3: revision and verification report',0)
report.add_paragraph('Main output: Astera_manuscript_draft3_rarefied.docx. Companion outputs: supplementary addendum, analysis workbook, high-resolution/vector figures, and revision scripts. Draft 2 and earlier manuscript files were preserved.')
report.add_heading('Substantive corrections',1)
for t in [
'Date-inclusive models explain adjusted fractions of 14.0%, 13.1%, and 9.1% for 16S, ITS, and 18S. Unique management is non-identifiable because all till observations are on 28 November and all NT observations are on other dates.',
'The analytical cohort spans 21 November–1 December 2022, with only 0–10 and 10–30 cm depths. Broader collection years/depths are distinguished from this cohort.',
'Observed ASV totals are 9,526, 155, and 790. Earlier values (33,112, 162, 3,142) counted all source-object columns, including zero-only columns in this subset.',
'Alpha diversity now uses sample means across 200 without-replacement rarefactions: 2,802 reads (16S), 10,488 (ITS), and 1,829 (18S). All samples are retained, and BH adjusts six tests. Main Figure 9 adds exact rarefaction curves. A higher-depth sensitivity retains at least 95% per marker; all main and secondary BH q values exceed 0.05. CLR, relative abundance, Hellinger, PCA, and PCoA are now distinguished.',
'All nine main figures were rebuilt, including three-marker ordination plus cross-marker PCA; new dbRDA and variance partitioning displays were integrated. The ITS carbon network is referenced in the archival supplement as Figure S31.',
'Added missing Cutadapt and DADA2 references, date/ecological methods references, and corrected the malformed Jansa author citation.'
]:report.add_paragraph(t,style='List Bullet')
report.add_heading('Verification completed',1)
for t in checks+['All sample-level read totals, observed richness, and Shannon values recomputed from count exports and matched to audit tables.','Cross-marker feature selection, CLR transformation, standardization, and PCA reconstructed; PERMANOVA R² matched the saved result.','All 61 comment IDs have an individual response in the companion CSV and below.','Main tables rebuilt from saved full-precision CSVs; provenance retained. Figures and embedded-image dimensions checked.']:
 report.add_paragraph(t,style='List Bullet')
report.add_heading('Unresolved source-dependent items',1)
for t in [
'Draft 2 supplies field/laboratory protocol details (farm sampling, extraction kit, primers, platform, Cutadapt and DADA2). These were retained as author-provided information, not independently verified against laboratory records. DADA2 version, filtering/trimming thresholds, taxonomy reference database/version, classifier settings, and PCR details remain unprovided.',
'BioProject PRJNA1371968 was supplied in draft 2, but public retrieval could not be verified. Confirm that it contains these samples and the release status before submission.',
'Pre-QC and intermediate sequencing-stage read counts are not in the available inputs. The new workbook explicitly labels those fields unavailable; retained matrix totals are provided.',
'Carbon method, fraction, and units are unknown; C/N/pH disagree for 20 matched AMF samples. Source chemistry reconciliation remains necessary. No favored replacement values were selected.',
'Workbook profile grouping enables restricted exploratory permutations but does not establish independent farm replication. The full field hierarchy and exchangeability assumptions need author review. The NumPy/SciPy dbRDA was rerun and algebraically validated, but not cross-validated in R/vegan.',
'The older main PERMANOVA, MaAsLin2, and network inference outputs were checked against saved tables; they were not all refitted from raw reads. All earlier archival supplementary images remain unchanged, with revised main displays superseding corresponding old figures.',
'The sensitivity figure can now be reviewed clearly, but the requested lab discussion (comment 59) remains a coauthor activity.',
'Word document structure and figure assets were checked programmatically. A native Microsoft Word/LibreOffice pagination renderer is unavailable in this environment; final pagination should be reviewed in Word. No guarantee of a mistake-free or submission-ready manuscript is made.'
]:report.add_paragraph(t,style='List Bullet')
report.add_heading('Comment-by-comment responses',1)
for _,r in a.iterrows():
 report.add_heading(f'Comment {r.comment_id}: {r.status}',2);report.add_paragraph(r.comment);report.add_paragraph(r.response)
report.save(O/'Draft3_revision_report.docx')
(O/'README.md').write_text('''# Draft 3 manuscript revision

Open `Astera_manuscript_draft3_rarefied.docx` for the revised manuscript. Read `Draft3_revision_report.docx` for all 61 comment responses, numerical checks, and unresolved source-verification items. The addendum extends the earlier supplement with Figure S54 and Tables S13–S15. `draft3_analysis_tables.xlsx` contains new analyses and sample-level retained read totals. All revised figure assets are in `figures/` (PNG, PDF, SVG).

Earlier manuscripts and outputs were not overwritten. The old supplementary atlas remains archival; use revised main figures for the updated displays. Do not call the draft submission-ready until the source-dependent items in the report are resolved.

## Reproduction
Run from the workspace root, using Python with numpy, scipy, pandas, statsmodels, matplotlib, threadpoolctl, openpyxl, python-docx, and lxml:

1. `python manuscript/draft3_rarefied/source/reproduce_dbrda.py`
2. `python manuscript/draft3_rarefied/source/rebuild_figures.py`
3. `python manuscript/draft3_rarefied/source/build_draft3.py`
4. `python manuscript/draft3_rarefied/source/review_report.py`

Step 2 also runs rarefied_alpha.py to regenerate alpha figures, exact curves, and all 200-draw sensitivity outputs. It can be run separately before rebuilding the Word files. The scripts read the existing source count exports and model outputs. The revision adds no raw-read reprocessing and does not claim R/vegan cross-validation.
''')
print('All 61 comments mapped; dbRDA results reproduced; revision report written.')
