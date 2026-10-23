# Astera computational package

Companion to *Tillage impacts on soil microbiomes differ by organismal group, community metric, and sensitivity to sampling imbalance.* Prepared 2026-10-23; not yet uploaded as a corrected public release.

## What changed
The companion corrected computational package adds missing final-revision outputs, shared model metadata with pseudonymous profile IDs, portable analysis scripts, numerical validation, and a final-result workflow index. Existing raw FASTQ archives are unchanged. Recorded study inputs remain 44 16S, 43 ITS and 44 AMF libraries, comprising 45 distinct sample IDs and 42 shared IDs. No metagenomic data are included.

## Run the validated workflows
From this extracted folder, using Python 3.9 and the recorded package versions:

```bash
python -m pip install -r environment/requirements-python-validated.txt
python source/run_final_checks.py
```

This runs processed-count checks, rarefied alpha diversity and unique-ASV retention, full-model dbRDA/four-set partitioning, three-set partitioning, profile-restricted PERMANOVA/date models, recorded/exclusion adjusted-carbon analyses, and manuscript figure generation. It compares the numerical outputs with archived reference tables. Outputs go to reproduced/ and validation/; final_reference/ stays unchanged.

The figure builder uses archived original management-only MaAsLin2, unrestricted PERMANOVA, balance-sensitivity and original network results for some panels. It does not refit those R models. source/prepare_R_inputs.R is an unvalidated adapter for reconstructing legacy phyloseq inputs; it is not part of the validated command. Original R scripts remain under source/original. Raw FASTQ-to-ASV processing and all historical R models have NOT been validated end-to-end by this release.

## Layout
- data/: selected counts, sequences, taxonomy, metadata, and METADATA_NOTES.md.
- final_reference/: final-revision results and figures, with source profile labels replaced by pseudonyms.
- reproduced/: regenerated results from the portable workflow (generated locally; not necessary as duplicated upload contents).
- results/ and figures/: original/legacy outputs retained for provenance.
- source/final/: portable final-analysis code; source/original/: archived original scripts.
- legacy_documentation/: original supplementary workbook and figure index, explicitly historical.
- FINAL_WORKFLOW_INDEX.csv: final outputs and generating commands.
- validation/: numerical comparison report and execution record.
- environment/: package versions and dependency guidance.

## Scientific status
Numerical reproduction is not validation of the sampling design or causal interpretation. Management and date remain confounded. C_16S/C_ITS/C_AMF retain discordant recorded chemistry; units/fraction and joins need author reconciliation. profile_id encodes existing grouping and does not establish field independence. The manuscript's adjusted-MaAsLin2 and date-restricted-PERMANOVA wording must be reconciled with actual analyses. Corrected manuscript text is not included here.

## Related archives and code
The separate denosing.zip and additional_scripts.zip hold archived processing and historical analysis scripts; their original paths and incomplete commands are documented. The code repository is https://github.com/Abdelrahim-maker/astera-microbiome-scripts . The earlier Zenodo version is https://doi.org/10.5281/zenodo.23582246 ; this local correction does not change that published record. Confirm licenses and author metadata before a new release.

For GitHub use, extract the corrected Zenodo processed ZIP into a working folder and place this repository content alongside its data/ and final_reference/ folders. The data and reference outputs are not stored in GitHub.
