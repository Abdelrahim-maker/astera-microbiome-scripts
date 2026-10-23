# Focused Astera microbiome release — 2026-10-23

Companion to *Soil microbial diversity and association networks under tillage and no-tillage in Ohio soybean fields*.

Retained: rarefied alpha diversity; marker-specific composition; management-only profile-restricted PERMANOVA; original exploratory management-only MaAsLin2; simple marker-specific network summaries; 20 balanced-subsample checks. No carbon, cross-marker, date-effect, dbRDA or variance-partition analyses are included. Collection dates remain as metadata documenting management/date confounding.

## Reproduce
Extract Astera_focused_data_results.zip and Astera_focused_code.zip into the SAME directory (both archives have relative paths, without a wrapper folder). Use Python 3.9 and install `environment/requirements-python-validated.txt`, then run:

```sh
python source/run_retained_checks.py
```

The runner recomputes rarefaction, the management-only restricted PERMANOVA, and six marker-specific networks from deposited counts; checks agreement with retained results; and rebuilds Word documents, figures and the workbook. It uses archived MaAsLin2, unrestricted/PERMDISP and balanced-subsample outputs when assembling documents. It does not refit those R analyses. A clean installation of all pinned dependencies has not been tested. `validation/` records executed checks.

Optional R adapters: `source/run_maaslin2_management.R` and `source/run_balanced_sensitivity.R`. They were not executed here because R/Maaslin2 is unavailable. They use only management, never date. The balance adapter records sample identities for NEW draws; the original subset-ID manifests were not recovered, and changing the historical workflow changes RNG consumption. Do not substitute its results for the archived tables without updating all related text and validating the rerun. Install Maaslin2 from Bioconductor and vegan/Hmisc/igraph/tibble through their documented package channels. No historical R lockfile is claimed.

## Files
`data/`: processed counts, ASV FASTA, taxonomy and metadata (no chemistry).
`inputs/`: retained numerical outputs and figure sources used by the document builder.
`manuscript/`: matching documents, figures, editable source and CSV tables.
`source/`: retained Python workflows and explicitly unvalidated R rerun adapters.
`provenance/model_logs/`: original marker-specific MaAsLin2 logs.
`validation/`: validation reports, with limits stated explicitly.
Raw reads are separate Astera_16S_raw_FASTQ.tar, Astera_ITS_raw_FASTQ.tar and Astera_AMF_raw_FASTQ.tar files on Zenodo. `denosing.zip` is historical processing provenance, not an end-to-end validated pipeline.

## Interpretation and provenance
Marker counts: 44/43/44; 45 distinct sample IDs, 42 shared. Profile IDs are pseudonymous and do not prove farm-level independence. Management is confounded with sampling date. PERMANOVA exchanges profiles within site and identical depth coverage; those restrictions are not fitted covariates. MaAsLin2 is unadjusted and exploratory. Network edges are reconstructed associations, not established interactions; agreement of graph totals does not prove identity with unavailable original edge exports. Historical raw-processing lineage and database versions remain incompletely documented.

## Release status and licenses
Prepared locally; no upload, new DOI or Git tag is claimed. Existing records: https://doi.org/10.5281/zenodo.23582246 and https://github.com/Abdelrahim-maker/astera-microbiome-scripts . Replace links with the new version DOI and immutable code tag after publication.

Proposed publication licenses: MIT for project code, CC0-1.0 for data/results, CC-BY-4.0 for manuscript text and figures; preserve third-party notices. These are release recommendations, not a claim that all coauthor rights have already been confirmed. Confirm the license declarations when publishing. See https://astera.org/open-science-policy/ .
