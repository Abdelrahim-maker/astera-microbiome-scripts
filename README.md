# Astera soil microbiome analysis scripts

Code accompanying **Tillage impacts on soil microbiomes differ by organismal group, community metric, and sensitivity to sampling imbalance.**

The manuscript cohort contains 44 16S, 43 ITS and 44 AMF samples. This repository contains all analysis and validation scripts from the prepared Zenodo resource, plus environment instructions and workflow documentation. Raw reads and research data belong in the companion Zenodo deposit; its DOI will be added after publication.

## Contents
- `source/original/`: original R analysis scripts and Python ordination scripts.
- `source/check_processed_data.py`: verifies sample counts and retained read totals.
- `source/reconstruct_network_edges.py`: reconstructs eight network edge tables and checks archived degrees/hubs.
- `environment/`: validated Python package pins and an R dependency installation recipe.
- `WORKFLOWS.md` and `WORKFLOW_INDEX.csv`: workflow guidance and historical provenance gaps.

## Use with the Zenodo data
Extract `Astera_data_sequences_scripts.zip` from the companion deposit into this repository root so that `data/` and `results/` sit alongside `source/`. Then run:

```bash
python -m pip install -r environment/requirements-python-validated.txt
python source/check_processed_data.py
python source/reconstruct_network_edges.py
```

Original R and ordination scripts retain their original workspace paths and need additional inputs described in WORKFLOWS.md. A fully self-contained rerun of every historical analysis is not claimed. The original R package-version lock and some historical postprocessing scripts were not found. Python pins describe the environment used for reconstruction, not the original R environment.

## Data and publication links
The GitHub URL and Zenodo DOI will be cross-linked after publication. No DOI is assigned by this repository.

## License
A software license has not yet been selected by the authors. No additional reuse license is granted by this preparation step.

## Funding
Astera Institute, Hidden Science Competition Award (2026).

## Additional original scripts
`denosing/` contains original amplicon-processing scripts; `additional_scripts/` contains historical manuscript and analysis revisions. See their READMEs and manifests for scope and limitations. These additions are archival originals, not a validated end-to-end pipeline.
