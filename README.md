# Astera soil microbiome analysis scripts

Code accompanying *Tillage impacts on soil microbiomes differ by organismal group, community metric, and sensitivity to sampling imbalance.*

This repository provides analysis scripts and workflow documentation for comparisons of no-till and tillage management across prokaryotic (16S), fungal (ITS), and arbuscular mycorrhizal fungal (AMF) communities. The manuscript cohort comprises **44 16S, 43 ITS, and 44 AMF samples**.

## Data and code

- **Zenodo data and computational materials:** https://doi.org/10.5281/zenodo.23582246
- **GitHub code:** https://github.com/Abdelrahim-maker/astera-microbiome-scripts

The Zenodo deposit includes paired-end raw reads, processed ASV abundance tables, representative sequences, taxonomy, sample metadata, statistical results, figures, and computational materials. Raw FASTQ reads are provided in three marker-specific TAR archives. No metagenomic materials are included.

## Repository contents

- `source/original/`: original R analysis scripts and Python ordination scripts.
- `source/check_processed_data.py`: verifies cohort sizes and retained read totals.
- `source/reconstruct_network_edges.py`: reconstructs eight network edge tables and checks archived node degrees and hub classifications.
- `denosing/`: original DADA2 and related amplicon-processing scripts for 16S, ITS and AMF, including QIIME 2 commands and R scripts.
- `additional_scripts/`: further manuscript, figure-generation, audit, dbRDA, carbon-analysis and rarefaction scripts, including historical revisions.
- `environment/`: Python package-version records and R dependency-installation instructions.
- `WORKFLOWS.md` and `WORKFLOW_INDEX.csv`: workflow guidance, result-to-script mappings and provenance notes.
- `SOURCE_MANIFEST.csv`: file inventory with SHA-256 checksums.

## Use with the Zenodo data

Download the processed-data archive from [Zenodo](https://doi.org/10.5281/zenodo.23582246). It may be named `Astera_data_sequences_scripts.zip` or `Astera_data_sequences_scripts (1).zip`. Extract it into a separate working folder, then place the repository contents in the same folder so that `data/`, `results/`, `source/` and `environment/` share the same root. Keep a copy of the original downloaded archives.

From that working folder, run:

```bash
python -m pip install -r environment/requirements-python-validated.txt
python source/check_processed_data.py
python source/reconstruct_network_edges.py
```

These commands check processed sample counts and read totals, then regenerate network edges and compare node degrees and hub classifications against archived results. They do not run all manuscript analyses.

For raw-read processing, consult `denosing/README.md`. The original scripts retain OSC paths, account settings and software setup commands that must be adapted to the execution environment. The ITS workflow uses primer-trimmed inputs. Historical script variants should not be executed sequentially without checking their intended inputs and cohort.

## Reproduction status

The data and scripts are publicly available to support reuse and independent evaluation of the manuscript analyses. Packaged count totals and reconstructed network node degrees and hub classifications have been checked against the archived outputs.

A clean end-to-end rerun of all final manuscript analyses using only the deposited materials has **not yet been verified**. Some scripts retain original workspace dependencies, and the exact original R package-version lock is unavailable. The Python package pins describe the environment used for reconstruction and validation. See `WORKFLOWS.md`, `WORKFLOW_INDEX.csv` and the folder READMEs for specific requirements and provenance gaps.

Saved PCoA coordinates and original NMDS figures are included in the Zenodo materials; they represent different ordination methods. Reconstructed network edge tables are identified as regenerated outputs. Historical and alternative-cohort analyses in `additional_scripts/` should be distinguished from the final manuscript analyses.

## Study limitations

Sampling date and management are confounded, limiting causal attribution of community differences to tillage. Carbon and other source-metadata discrepancies are documented in the accompanying materials. Co-occurrence networks represent exploratory associations.

## License

Consult the repository's `LICENSE` file, when present, for software reuse terms. Public availability alone does not grant a software license. Dataset reuse terms are specified on the [Zenodo record](https://doi.org/10.5281/zenodo.23582246); software and dataset licenses should be considered separately.

## Funding

Astera Institute, Hidden Science Competition Award (2026).
