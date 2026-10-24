# Astera soil microbiome analysis

Sequencing source data and analysis code supporting:
“Soil microbial diversity and association networks under tillage and no-tillage in Ohio soybean fields.”

## Current focused package

Start with the [focused package README](github_upload/README.md).

- [Analysis scripts](github_upload/source/)
- [Computational environment](github_upload/environment/)

The focused release contains sequencing source data and analysis code. It excludes manuscript text, supplementary documents, figures, workbooks and archived manuscript result tables.

The manuscript and supplementary document are maintained separately for bioRxiv.

## Data release

Published Zenodo dataset, version 2.0.0:
https://doi.org/10.5281/zenodo.23603407

Earlier dataset version:
https://doi.org/10.5281/zenodo.23582246

Extract `Astera_source_data.zip` into the same working directory as the contents of `github_upload/`. Follow the focused package README for execution instructions.

## Analysis scope and reproducibility

The Python workflow calculates rarefaction, management-only profile-restricted PERMANOVA and marker-specific networks from source inputs. Running the code generates new numerical outputs locally; no figures are generated.

Optional management-only MaAsLin2 and balanced-subsampling R adapters are provided but were not executed in the preparation environment. The balanced adapter generates new subsets; it does not recover the historical selections.

Collection dates remain metadata, not model predictors. Carbon, cross-marker and date-effect analyses are excluded from the focused workflow.

Management and collection date are confounded. Recorded profile identifiers do not establish independent farm replication. The complete historical raw-read-to-final-table processing lineage remains incompletely documented.

## Historical materials

[astera-microbiome-scripts-github-corrected/](astera-microbiome-scripts-github-corrected/) preserves the earlier computational package and is superseded for the current focused release.

Root-level `WORKFLOWS.md`, `WORKFLOW_INDEX.csv` and `SOURCE_MANIFEST.csv` describe earlier materials. Use the focused package README for the current workflow.

## Licensing

Project-authored code is licensed under the [MIT License](LICENSE). Dataset licensing is specified separately in the Zenodo record.
