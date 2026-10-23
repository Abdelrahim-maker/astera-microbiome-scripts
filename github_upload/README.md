# Astera sequencing data and analysis code

This focused release contains sequencing source data and analysis code only. It excludes the full manuscript, supplementary documents, workbooks, figures, captions, manuscript result tables, archived fitted-model results and network edge exports. The manuscript and supplement remain separate for bioRxiv.

## Data
Reserved Zenodo DOI: https://doi.org/10.5281/zenodo.23603407 (activates on publication).
Earlier dataset: https://doi.org/10.5281/zenodo.23582246 .
The source-data ZIP contains ASV counts, representative FASTA sequences, taxonomy, and sample metadata. Three raw FASTQ archives and historical denosing.zip are separate files. These underlying sequencing measurements remain included; no derived manuscript results are deposited in this version.

## Run
Extract Astera_focused_code.zip and Astera_source_data.zip into the same directory. If using GitHub, use the contents of github_upload/ as this directory. With Python 3.9 and the dependencies in environment/requirements-python-validated.txt:

```sh
python source/run_retained_checks.py
```

This runs rarefaction, management-only whole-profile PERMANOVA and marker-specific network reconstruction. New numerical outputs are written locally under reproduced/. No figures or preprint documents are generated. There is no comparison against archived manuscript tables because those tables are intentionally excluded.

Optional R adapters source/run_maaslin2_management.R and source/run_balanced_sensitivity.R are provided but were not executed in this environment. The balanced adapter creates new subsets and records their IDs; it does not recover historical subsets. MaAsLin2 is management-only. The Python runner does not refit either R analysis. No historical R lockfile or validated end-to-end raw-read pipeline is claimed.

## Scope and limitations
Code analyzes each marker separately. Carbon, cross-marker and date-effect analyses are excluded. Dates remain metadata only. Management and collection date are confounded; recorded profiles do not establish farm-level independence. Networks are statistical associations rather than demonstrated interactions. The full raw-read-to-final-table lineage remains incompletely documented.

## Licensing
Publication licenses remain to be finalized. LICENSE_PLAN.md is a plan, not a license grant. Preserve third-party notices.
