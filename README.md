# Astera soil microbiome analysis

Computational materials for:
“Soil microbial diversity and association networks under tillage
and no-tillage in Ohio soybean fields.”

## Current manuscript and workflows

The current focused package is in [github_upload/](github_upload/).

Start with its [README](github_upload/README.md).

- [Analysis scripts](github_upload/source/)
- [Computational environment](github_upload/environment/)
- [Manuscript sources](github_upload/manuscript/source/)

Retained analyses:
- Rarefied alpha diversity and marker-specific community composition.
- Profile-restricted, management-only PERMANOVA.
- Exploratory management-only MaAsLin2 associations.
- Simplified marker-specific networks.
- Balanced-subsampling sensitivity analyses.

Collection dates remain metadata, not model predictors.

## Historical materials

[astera-microbiome-scripts-github-corrected/](astera-microbiome-scripts-github-corrected/)
preserves the earlier package and is superseded for the revised manuscript.

The root-level WORKFLOWS.md, WORKFLOW_INDEX.csv and SOURCE_MANIFEST.csv
also describe earlier materials. Use the focused package README
for the current workflow.

## Data release

 Zenodo version: https://zenodo.org/records/23603407

The matching focused Zenodo version is pending. The earlier deposit
should not be assumed to contain the revised package.

## Reproducibility

The focused workflow recomputes rarefaction, restricted management
PERMANOVA and marker-specific network summaries. It uses archived
MaAsLin2 and historical balanced-subsample results; those R analyses
were not refitted.

Management and collection date are confounded. The complete
raw-read-to-final-table processing lineage remains incompletely documented.

## Licensing

Project-authored code is licensed under the [MIT License](LICENSE).
Third-party materials retain their existing licenses.
Dataset licensing is specified separately in the Zenodo record.
