# Astera focused analysis package

This directory contains the public, GitHub-ready focused release for the Astera microbiome analysis.

The package is intended to make the workflow easy to inspect, rerun, and audit while keeping the publication documents and manuscript artifacts separate from the code and data release.

## Package contents

- [source](source) — focused analysis workflow and validation scripts
- [environment](environment) — validated environment metadata
- [manuscript](manuscript) — manuscript generation and validation support
- [reproducibility_update_2026_10_27](reproducibility_update_2026_10_27) — public reproducibility package and audit results

## Public data release

This repository references the public release records:

- Zenodo v2.0.0: https://doi.org/10.5281/zenodo.23603407
- Earlier release: https://doi.org/10.5281/zenodo.23582246

The source-data archive and focused analysis code are intended for transparent reuse. The manuscript and supplementary documents are maintained separately from the code repository.

## Quick start

Use the validated Python environment described in [environment/requirements-python-validated.txt](environment/requirements-python-validated.txt).

From this folder, run:

```bash
python source/run_retained_checks.py
```

This script executes the retained analysis workflow for the focused package, including rarefaction, management-only PERMANOVA, and marker-specific network reconstruction.

## Reproducibility and audit package

The public reproducibility update is here:

- [reproducibility_update_2026_10_27](reproducibility_update_2026_10_27)

This folder contains:

- the packaged archive
- verification notes
- comparison metadata
- execution metadata
- balanced comparison output
- R environment details

This material is intended to make the audit trail easy to inspect and share publicly without overstating the scope of verification.

## Scope and limitations

The focused package is intentionally scoped to the retained analyses and does not claim a full public reconstruction of every historical workflow step.

Important limitations include:

- MaAsLin2 remains dependent on the required package installation.
- The balanced adapter creates new subsets and does not recover the historical selection.
- Management and collection date remain confounded in the study design.
- The complete raw-read-to-final-table lineage remains incompletely documented.

## Licensing

This repository is distributed under the [MIT License](../LICENSE), while dataset licensing remains as specified in the associated Zenodo record.
