# Soil microbial diversity and association networks under tillage and no-tillage in Ohio soybean fields

Afaf Abdelrahim¹, Alec Ogg³, Klaus Lorenz², Rattan Lal², M. Scott Demyan³, Virginia I. Rich¹

¹ Department of Microbiology, The Ohio State University, Columbus, OH, USA

² CFAES Rattan Lal Center for Carbon Management and Sequestration, The Ohio State University, Columbus, OH, USA

³ Department of Horticulture and Crop Science, The Ohio State University, Columbus, OH, USA

Corresponding author: Virginia I. Rich; rich.270@osu.edu

## Abstract

Tillage can influence soil microbial communities, but observed associations may differ among organismal groups and analytical measures. We compared prokaryotes, general fungi, and arbuscular mycorrhizal fungi (AMF) in soils from commercial soybean fields in Ohio under no-tillage (NT) and tillage. The analytical datasets comprised 44 16S rRNA gene libraries, 43 fungal ITS libraries, and 44 AMF-targeted 18S rRNA gene libraries from two soil depths. Rarefied richness and Shannon diversity did not differ between management groups after false-discovery-rate correction. Management explained 4.7%, 6.5%, and 3.1% of Bray–Curtis variation in 16S, ITS, and AMF communities, respectively; none of the profile-restricted PERMANOVA tests met the adjusted significance threshold (q = .0555, .3770, and .3018). Exploratory management-only MaAsLin2 models identified 5, 14, and 1 associated ASVs, respectively, at q ≤ .05. Correlation-network topology differed between groups, but several contrasts changed substantially after balancing sample numbers, particularly AMF connectivity. These results identify candidate taxon associations and illustrate the sensitivity of community comparisons to sampling structure. Because management was confounded with collection date and the study was observational, the findings do not establish independent or causal effects of tillage.

Keywords: soil microbiome; no-tillage; tillage; alpha diversity; beta diversity; arbuscular mycorrhizal fungi; differential abundance; association networks

## Introduction

Soil microorganisms contribute to nutrient cycling, organic matter decomposition, and plant–soil interactions [@Fierer2017]. Agricultural management can change their habitats, including the distribution of plant residues and the physical continuity of soil pores and fungal hyphae. Community responses may vary with soil depth and among microbial groups [@Sun2018]. Studies of cropping practices have documented changes in bacterial and fungal community composition [@Hartman2018; @Sommermann2018], while a meta-analysis identified management-related responses of arbuscular mycorrhizal fungi (AMF) [@Bowles2017]. These findings motivate comparisons that include multiple microbial groups in working agricultural fields.

Different analytical measures address different questions. Alpha diversity summarizes variation within samples, beta diversity describes compositional differences among samples, and differential-abundance models identify individual features associated with recorded management. Correlation networks summarize statistical associations among features, but their topology is sensitive to filtering, sample size, and the compositional nature of sequencing data [@Gloor2017; @Rottjers2018; @Weiss2016]. Consequently, agreement across these measures should not be assumed.

We compared no-tillage (NT) and tillage soils from commercial soybean fields in Coshocton County, Ohio, using 16S rRNA gene, fungal ITS, and AMF-targeted 18S rRNA gene amplicons. We examined alpha diversity, community composition, management-associated ASVs, and association networks. We retained checks for sampling structure and unequal management-group sizes to evaluate the robustness of these comparisons. Our objective was to describe management-associated patterns in the sampled soils, rather than infer causal tillage effects.

## Materials and methods

### Soil sampling and analytical cohorts

Soils were collected from commercial soybean-growing farms within six watersheds in Coshocton County, Ohio, on predominantly silt loam soils. Management histories were reported as long-term no-tillage, maintained for at least 40 years, or tillage to approximately 20 cm. Sampling polygons represented uniform soil series and management histories. Within polygons, three sampling centroids were selected at least 20 m apart. Soil was collected at each centroid and at four surrounding points approximately 10 m away using a 3.2-cm-diameter battery-powered auger, then composited by depth interval. The retained depth intervals were 0–10 and 10–30 cm. Samples were transported in coolers to the laboratory within 8 h and stored at −80 °C.

Collection occurred on November 21, 22, and 28 and December 1, 2022. All tillage samples were collected on November 28, and no NT samples were collected that day. Management and collection date therefore cannot be separated in these data. Farmer identities and exact locations are withheld to protect privacy. Pseudonymous profile identifiers preserve the recorded grouping of samples across depths; the available grouping records do not establish independent farm-level replication.

The retained matrices contained 44 16S libraries (33 NT, 11 tillage), 43 ITS libraries (33 NT, 10 tillage), and 44 AMF libraries (34 NT, 10 tillage). These correspond to 45 distinct sample identifiers across markers, with 42 shared identifiers (33 NT, 9 tillage), rather than one identical 44-sample cohort. Marker membership, depth, recorded profile grouping, site code, and collection date are provided in Table S1. Absence from a retained matrix is not treated as evidence of a particular sequencing failure.

### DNA extraction, sequencing, and processed inputs

DNA was extracted from 0.25 g of soil using the DNeasy PowerSoil Pro Kit (Qiagen). Amplicon preparation and Illumina MiSeq sequencing were performed by the Environmental Sample Preparation and Sequencing Facility at Argonne National Laboratory. The reported primer pairs were 515F/806R for prokaryotic 16S rRNA genes [@Caporaso2012], ITS1F/ITS2 for general fungi [@Gardes1993; @White1990], and NS31/AML2 for AMF-targeted 18S rRNA genes [@Lee2008; @Simon1992]. The ITS and AMF-targeted assays were analyzed separately and were not treated as independent estimates of the same fungal community.

Archived processing workflows include primer trimming with Cutadapt and amplicon denoising with DADA2 through run-specific R or QIIME2 implementations [@Callahan2016; @Martin2011]. Analyses used the retained ASV count tables, representative sequences, and associated taxonomic annotations organized in phyloseq objects [@McMurdie2013]. The documented reference databases were SILVA for 16S, UNITE for ITS, and MaarjAM for AMF [@Nilsson2019; @Opik2010; @Quast2013]. Exact run-to-final-object provenance and database releases are incompletely documented; the retained processing records and these limitations are summarized in Supplementary Methods. No uniform trimming configuration or complete raw-read-to-result reproduction is assumed.

### Alpha diversity

Counts were rarefied without replacement to the minimum retained library size for each marker: 2,802 reads for 16S, 10,488 for ITS, and 1,829 for AMF. Richness and Shannon entropy, using natural logarithms, were averaged over 200 draws per sample. NT and tillage sample means were compared using two-sided asymptotic Mann–Whitney U tests, equivalent to Wilcoxon rank-sum tests, with tie and continuity corrections. Benjamini–Hochberg (BH) adjustment was applied across the six marker-by-metric comparisons [@Benjamini1995]. Rarefaction draws were not treated as independent observations. These sample-level tests do not adjust for profile dependence or design covariates.

A sensitivity analysis used higher depths retaining at least 95% of samples: 3,799, 16,867, and 5,882 reads, respectively. Exact expected-richness curves describe subsampling without extrapolation [@Hurlbert1971]. Sample means, sensitivity tests, and library depths are provided in Table S2 and Supplementary Data; curves are shown in Figure S1.

### Beta diversity and management tests

Bray–Curtis dissimilarities were calculated from marker-specific relative abundances. Principal coordinates analysis (PCoA) with Lingoes correction was used for visualization; management tests used the original Bray–Curtis distances. PERMANOVA tested the model Bray–Curtis distance ~ management [@Anderson2001]. Site code and depth coverage defined permutation blocks; they were not fitted covariates. Collection date was neither a predictor nor a permutation stratum. For the main comparison, 9,999 permutations exchanged whole recorded soil profiles within the same site code and identical depth coverage, preserving pairing and depth positions. BH adjustment was applied across the three marker tests. This scheme accommodates the recorded profile grouping but does not isolate management from collection date or establish farm-level independence.

Original management-only PERMANOVA tests with 999 unrestricted sample-level permutations are reported separately as exploratory comparisons in Table S3. Dispersion was assessed with the archived PERMDISP analysis using 999 permutations [@Anderson2006]. A nonsignificant dispersion test does not establish equal dispersions, particularly with unbalanced groups [@Anderson2013]. To assess sensitivity to unequal group sizes, 20 overlapping subsets retained all tillage samples and randomly selected the same number of NT samples: 11 per group for 16S and 10 per group for ITS and AMF. These original sensitivity tests used unrestricted permutations and were not profile-matched; detection frequencies are descriptive, not independent replications or power estimates.

### Differential abundance

Management associations were evaluated separately for each marker using MaAsLin2 linear models [@Mallick2021]. ASVs present in at least 10% of samples, with a minimum of three samples, were retained before dividing counts by each sample’s total over the retained features to obtain relative abundance. The recorded MaAsLin2 settings were normalization = NONE, transform = LOG, analysis_method = LM, and correction = BH. The model formula was transformed relative abundance ~ management, with NT as the reference. No site, depth, collection-date, or random-profile effects were included. These are unadjusted exploratory models, not covariate-adjusted estimates of tillage effects. Positive coefficients therefore indicate higher modeled transformed relative abundance under tillage.

BH correction was performed within each marker. We report q ≤ .05 associations and separately summarize q ≤ .25 exploratory candidates. These thresholds refer to the fitted management-only models, whose assumptions do not account for the field hierarchy. Complete results, including nonsignificant tests, are supplied in Supplementary Data.

### Association networks

Networks were inferred separately for each marker and management group. Candidate ASVs had prevalence ≥10% within the group, and at most 250 were retained by summed abundance. The source workflow used the 100 most abundant features if fewer than 10 passed the prevalence filter. A pseudocount of one was added before centered log-ratio transformation [@Gloor2017]. Spearman correlations were retained as edges when |ρ| ≥ .6 and BH-adjusted q < .05 across tested feature pairs. Only ASVs incident to retained edges were included as graph nodes.

We summarized only connected ASV counts, edge counts, density, and mean degree. Density is the fraction of possible ASV pairs connected by an edge, and mean degree is the mean number of connections per connected ASV. These summaries avoid assigning ecological roles to individual nodes. Networks were reconstructed for the 20 balanced subsets described above. Network summaries are descriptive: graph nodes and resampled subsets are not independent field replicates, and correlations do not demonstrate direct interactions or ecological function [@Rottjers2018; @Weiss2016].

## Results

### Analytical coverage and alpha diversity

The retained cohorts contained 9,526 observed ASVs for 16S, 155 for ITS, and 790 for AMF, with median retained read counts of 7,611, 28,678, and 15,417.5, respectively (Table 1). These are processed-table totals, not raw sequencing yields. Rarefaction curves describe diversity within the retained matrices and do not demonstrate complete recovery of soil diversity (Figure S1).

Neither rarefied richness nor Shannon diversity differed between NT and tillage after BH correction (q ≥ .322; Figure 1; Table 2). AMF richness and Shannon diversity were numerically higher in tillage samples, but neither met the adjusted threshold. The higher-depth sensitivity analysis reached the same overall conclusion (Table S2).

### Community composition depended on the permutation scheme

Management accounted for 4.7%, 6.5%, and 3.1% of Bray–Curtis variation in the 16S, ITS, and AMF datasets, respectively (Figure 2; Table 3). Whole-profile permutations restricted by site code and depth coverage yielded p = .0185, .3770, and .2012, respectively; after BH adjustment, none met q ≤ .05 (q = .0555, .3770, and .3018).

The original unrestricted tests yielded p = .006 for both 16S and ITS and p = .072 for AMF (Table S3). No dispersion difference was detected at p < .05 (16S, p = .565; ITS, p = .481; AMF, p = .074). In balanced sample-level subsets, the unrestricted PERMANOVA p value was below .05 in 18 of 20 iterations for 16S, 7 for ITS, and 4 for AMF (Figure S2; Table S7). These sensitivity results do not supersede the restricted tests or resolve the sampling confounding.

### Management-associated ASVs

MaAsLin2 tested 915 16S, 127 ITS, and 87 AMF features. At q ≤ .05, 5, 14, and 1 ASVs, respectively, were associated with management; at q ≤ .25, the corresponding totals were 34, 43, and 3 (Figure 3). The broader totals include the q ≤ .05 features. Table S4 provides all 20 q ≤ .05 associations, and Supplementary Data includes all 1,129 tested features.

Among the five 16S associations, ASVs assigned to *Candidatus Xiphinematobacter*, Rokubacteriales, and Acidobacteriota Subgroup 5 had negative coefficients, indicating higher modeled abundance in NT. Gemmataceae and *Ammoniphilus* ASVs had positive coefficients. Nine of the 14 ITS associations were higher in NT and five in tillage. NT-associated annotations included *Neocosmospora rubicola* and *Clonostachys rosea*, whereas a *Tricellula* ASV was associated with tillage. The single AMF association, assigned to *Glomus* cf. *microaggregatum*, was higher in tillage. These are feature-level annotations and exploratory associations, not verified functional traits or independent tillage indicators.

### Network contrasts were sensitive to sample balance

Original network mean degree was higher in NT than tillage for 16S (9.79 versus 3.84), similar for ITS (5.46 versus 5.11), and lower in NT for AMF (6.69 versus 42.54; Figure 4; Table S5).

Balancing sample numbers changed the magnitude, and for some metrics the direction, of the contrasts (Table S6). Median NT mean degree became 7.18 for 16S, 3.56 for ITS, and 36.34 for AMF; tillage values remained 3.84, 5.11, and 42.54 because all tillage samples were retained. Thus, the large AMF connectivity contrast contracted markedly, and the ITS mean-degree ordering reversed. These changes demonstrate that the original topology contrasts cannot be interpreted independently of sample number and feature selection.

## Discussion

The study provides a focused comparison of microbial diversity, composition, feature-level associations, and network topology in working soybean fields. Rarefied alpha diversity did not differ detectably between management groups. This does not establish equivalence, given the smaller tillage cohorts, sample dependence, and the limited diversity captured by retained amplicon tables. Alpha diversity and community composition address distinct aspects of community organization, so a lack of alpha-diversity differences does not imply identical taxonomic composition.

Management explained a modest proportion of compositional variation, but the inferential conclusion depended on the permutation scheme. The 16S and ITS associations detected by unrestricted tests did not pass BH correction when profiles were exchanged within site code and depth coverage. The restricted 16S result was close to, but above, the .05 adjusted threshold and should not be described as significant. Earlier studies demonstrate that cropping practices can be associated with microbial composition [@Hartman2018; @Sommermann2018], but those findings do not establish the cause of the patterns observed here.

The MaAsLin2 results identify features for follow-up rather than validated tillage indicators. The fitted models contained management alone and did not account for profile dependence, site, depth, or collection date. Different ASVs assigned to a genus need not share ecological roles, and taxonomic annotations alone do not demonstrate beneficial plant effects, pathogenicity, or functional activity. Independent sampling and targeted validation are needed to evaluate the identified associations.

Network comparisons were especially sensitive to sample balance. The AMF NT network became much more connected when inferred from smaller NT subsets, reducing the original contrast with tillage. This behavior is consistent with the broader methodological sensitivity of microbial correlation networks [@Rottjers2018; @Weiss2016]. Greater connectivity or lower modularity should not be interpreted here as greater community stability, cooperation, or resilience.

Several design and data limitations constrain interpretation. Management was confounded with collection date, and available profile identifiers do not establish independent farm replication. Sample numbers were unequal, marker membership differed, and the historical preprocessing records do not completely establish run-to-final-table provenance. These limitations remain even when results are numerically reproducible from the processed inputs. Balanced sampling across management, collection occasions, and independent fields would strengthen future comparisons.

## Conclusions

Rarefied alpha diversity showed no detected management association, and profile-restricted beta-diversity tests did not pass correction across markers. Exploratory MaAsLin2 models identified management-associated ASVs, while inferred network contrasts changed substantially with sample balance. Together, the results describe microbial patterns in the sampled NT and tillage soils and emphasize the importance of sampling structure when interpreting community comparisons.

## Data and code availability

The earlier data archive is available at https://doi.org/10.5281/zenodo.23582246 and the code repository is https://github.com/Abdelrahim-maker/astera-microbiome-scripts. A focused revision package containing the processed inputs, retained outputs, manuscript sources, and analysis code has been prepared locally for a new release; the earlier archive is not asserted to contain this revision. The version-specific archive DOI and code release identifier will be added after deposition. The accompanying Supplementary Data workbook contains the retained numerical results, including nonsignificant differential-abundance tests. This availability statement does not imply that every historical analysis or the complete raw-read processing workflow has been independently reproduced.

## Acknowledgments

We acknowledge the Carbon Farming (C-FARM) initiative, led by Rattan Lal at The Ohio State University, for the research framework and support for soil sampling. We thank the participating farmers in Coshocton County for access to their fields, the Environmental Sample Preparation and Sequencing Facility at Argonne National Laboratory for sequencing services, and the Ohio State University teams involved in field sampling, soil processing, and laboratory work.

## Funding

This research was conducted as part of the Carbon Farming Alliance for Research and Management (C-FARM), led by Rattan Lal at The Ohio State University and supported by the Foundation for Food & Agriculture Research and collaborating funding partners. Additional support for manuscript preparation and research dissemination was provided by the Astera Institute through the Hidden Science initiative.

## Supplementary material

The companion Supplementary Material contains Supplementary Methods, Tables S1–S7, and Figures S1–S2. Supplementary Data provides machine-readable metadata and statistical results for the retained analyses, including all MaAsLin2 tests.
