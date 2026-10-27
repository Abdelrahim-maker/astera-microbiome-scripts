# Astera microbiome analysis

A public-facing repository for the Astera soil microbiome analysis workflow, including the focused GitHub release and a reproducibility update package.

This project supports the analysis of soil microbial diversity and association networks under tillage and no-tillage in Ohio soybean fields.

## Overview

The repository is organized to make the project easy to navigate:

- [github_upload/README.md](github_upload/README.md) — focused package overview and execution notes
- [github_upload/source](github_upload/source) — analysis scripts and validation workflow
- [github_upload/environment](github_upload/environment) — validated environment metadata
- [github_upload/reproducibility_update_2026_10_27](github_upload/reproducibility_update_2026_10_27) — public reproducibility package and audit artifacts
- [SOURCE_MANIFEST.csv](SOURCE_MANIFEST.csv) — source manifest
- [WORKFLOW_INDEX.csv](WORKFLOW_INDEX.csv) — workflow overview

## Current public package

The primary public-facing analysis package is in [github_upload](github_upload).

Quick links:

- [Focused package README](github_upload/README.md)
- [Analysis scripts](github_upload/source)
- [Environment files](github_upload/environment)
- [Reproducibility package](github_upload/reproducibility_update_2026_10_27)

## Reproducibility update

A GitHub-ready reproducibility update is included in [github_upload/reproducibility_update_2026_10_27](github_upload/reproducibility_update_2026_10_27). It contains:

- the packaged reproducibility archive
- verification notes
- comparison metadata
- execution metadata
- balanced comparison results
- the R environment log

This package is intended to make the evidence trail and the local audit accessible to the public in a clean, easy-to-navigate format. It is a GitHub-ready update and should not be interpreted as a replacement for the original public Zenodo release.

## Data release and context

The project references the published Zenodo release and earlier dataset version here:

- Zenodo v2.0.0: https://doi.org/10.5281/zenodo.23603407
- Earlier release: https://doi.org/10.5281/zenodo.23582246

The repository focuses on processing scripts, environment metadata, and source-data handling rather than manuscript text or derived figure outputs.

## Quick start

For the focused analysis workflow, start with:

- [github_upload/README.md](github_upload/README.md)
- [github_upload/source/run_retained_checks.py](github_upload/source/run_retained_checks.py)

The code is designed to run locally and generate numerical outputs without producing manuscript figures.

## Scope and limitations

This project is organized for transparency and reproducibility, but it is important to keep the limitations explicit:

- The Python workflow covers the retained analyses in the focused package.
- The optional MaAsLin2 and balanced-subsampling adapters are provided but not treated as complete historical reproduction.
- Management and collection date are confounded in the study design.
- The historical raw-read-to-final-table lineage remains incompletely documented and is not claimed as fully reconstructed.

## Licensing

Project-authored code is licensed under the [MIT License](LICENSE). Dataset licensing remains as specified in the published Zenodo record.

## Contact and maintenance

This repository is intended to provide a clear public record of the Astera workflow, audit artifacts, and analysis package for reproducibility and review.
