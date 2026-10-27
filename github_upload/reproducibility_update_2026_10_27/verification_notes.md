# Reproducibility audit — 27 October 2026

## Scope and evidence

This audit freshly extracted the local release copies of `Astera_focused_code.zip`, `Astera_source_data.zip`, and `denosing.zip`. Original archives were preserved. Their SHA-256 hashes and actual Python versions are in `verification/execution.json`. We could not independently retrieve the Zenodo records; equivalence to the publicly hosted bytes and published manuscript is not independently verified.

Comparison references are the retained CSV tables in the local 26 October manuscript revision, copied unchanged into `verification/reference/`. “Verified” below means agreement with those explicitly identified tables, not verification of every result in the published manuscript or of raw-read provenance. Numeric tolerance: relative 1e-10, absolute 1e-12; row order, columns and nonnumeric values must also agree. No scientific parameters or manuscript results were changed to obtain agreement.

## Checklist

| Analysis | Status | Evidence / manuscript target |
| --- | --- | --- |
| 16S, ITS, AMF rarefaction, richness and Shannon | Verified against local references | Table 2, both alpha scenarios, depth selection, sample alpha, per-draw summary and unique-ASV retention; 6 tables match |
| Bray–Curtis management PERMANOVA | Verified against local reference | Table 3; 9,999 whole-profile permutations restricted by site/depth coverage; q = 0.0555 / 0.3770 / 0.3018 for 16S / ITS / AMF |
| MaAsLin2 differential abundance | Not yet tested successfully | Adapter attempted, exit 1: package `Maaslin2` unavailable. Targets of 5 / 14 / 1 ASVs at q ≤ 0.05 remain unverified |
| Balanced subsampling | Runs but differs | 20 subsets per marker, 60 beta rows, 120 network rows; equal NT/till sample counts and unique sample IDs checked; 16 of 39 historical summary values match. No ASV detection-frequency output is implemented |
| Microbial association networks | Verified against local reference | S4: all 6 networks match connected nodes, edges, density and mean degree. Individual edge identity is not compared with a historical edge reference |
| Raw FASTQ → DADA2 → final tables/taxonomy | Historical workflow incomplete | Scripts inspected, not executed; see HISTORICAL_DADA2.md |

Only the three locally verified analysis categories are complete within this limited audit. End-to-end public-deposit reproduction remains incomplete.

## Inputs, software and exact commands

Unzip the source-data archive into the code directory. Required inputs are `data/{16S,ITS,AMF}_counts.csv.gz`, matching `*_samples.csv`, and `data/sample_metadata.csv`; all are included in the source-data archive. Counts are samples by ASVs. The Python metadata adapter maps `profile_id` to `Site_Name`. Keep marker-specific sample IDs and row names unchanged.

The run used Python 3.9.21, NumPy 1.23.5, pandas 2.3.3, SciPy 1.9.3, statsmodels 0.14.6 and threadpoolctl 3.6.0. See `verification/execution.json`. A fresh installation of the full historical requirements file was not tested.

Executed from a fresh extracted code/data root:

```sh
python source/run_retained_checks.py
```

Exit code 0. It executes `source/final/rarefied_alpha.py`, `source/final/rarefied_unique_asvs.py`, `source/final/restricted_management.py`, then `source/reconstruct_marker_networks.py`. Rarefaction uses 200 draws per sample/scenario; parameters and seeds remain those in the release scripts.

The added comparison command, run from the project workspace, was:

```sh
python biorxiv_submission_2026_10_23/focused_release/github_upload/source/compare_reproduction.py --run-root Astera_reproducibility/analysis_code --report Astera_reproducibility/comparison.json
```

Exit code 0; all eight comparisons pass. For a newly extracted corrected GitHub folder after running the analyses, the equivalent portable command is:

```sh
python source/compare_reproduction.py --report reproduced/comparison.json
```

R execution used OSC modules `gcc/12.3.0` and `R/4.4.0`:

```sh
module load gcc/12.3.0 R/4.4.0
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript source/run_balanced_sensitivity.R
Rscript source/run_maaslin2_management.R
```

Balanced command exit 0; MaAsLin2 command exit 1 (missing package). Other systems should use their own R installation instead of OSC module commands. R package versions and session information are recorded in `verification/R_environment.log`; the successful balanced run also writes `reproduced/balanced_sensitivity/sessionInfo.txt`. Required balanced packages: vegan, Hmisc, igraph and tibble; differential abundance requires Maaslin2. These are observed dependencies, not an original R lockfile.

## Outputs and comparison results

- `reproduced/rarefied/rarefaction/`: per-draw values, sample means, alpha tests, depths, exclusions, expected-richness curves and ASV retention summaries.
- `reproduced/restricted_permanova/model_statistics.csv`: per-marker R², pseudo-F, p, q and permutation settings.
- `reproduced/network_edges/`: six generated edge tables and `network_summary.csv`.
- `reproduced/balanced_sensitivity/`: selected sample IDs, beta statistics, network statistics and R session.
- `verification/comparison.json`: eight Python reference comparisons.
- `verification/balanced_comparison.csv`: all 39 observed versus historical balanced summary values.

Balanced median PERMANOVA p-values were 0.0150 / 0.0575 / 0.1675 (16S / ITS / AMF), versus historical 0.0320 / 0.0765 / 0.1205. Historical subset membership is unavailable; the deposited adapter generates new draws and uses unrestricted sample permutations. It does not reproduce the profile-restricted full-cohort analysis. The mismatch is documented, not corrected by changing seeds or results. Its seed is 20260904.

## Known limitations

1. MaAsLin2 remains blocked by the missing dependency. Even after installation, compare the actual q ≤ 0.05 results with all reference tests; the adapter's output threshold of 0.25 is not the manuscript discovery threshold.
2. The balanced adapter does not refit MaAsLin2 or produce per-ASV detection frequencies. Its successful beta/network run is not verification of the complete historical sensitivity workflow.
3. Original historical subsets and the R environment lock were not recovered. New balanced draws cannot establish exact historical reproduction.
4. Raw-read lineage, primer trimming and database checksums remain incomplete. Historical scripts contain conflicting variants and absolute paths; see HISTORICAL_DADA2.md.
5. Full clean dependency installation, historical edge identities, public-archive byte identity and direct published-manuscript comparisons were not verified.
6. Management/date confounding, profile dependence and provisional AMF taxonomy remain interpretation limitations.

## Publication handling

This is a local GitHub-ready update. No remote repository or Zenodo record was changed. The original Version 2.0.0 ZIPs remain untouched. README/reproducibility updates can be committed to GitHub. If the changed code archive or new reference tables are deposited on Zenodo, publish a new version and retain the original release. Do not describe these added files as already present in Version 2.0.0.
