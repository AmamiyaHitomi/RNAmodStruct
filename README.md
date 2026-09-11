# RNAmodStruct (ModStruct)

Human m6A modification level vs. local RNA structure signal: association and predictive-value analysis
across HEK293T (discovery) and HeLa (replication), using public GLORI (m6A) and icSHAPE (structure) data.

## Environment

- Python: `E:\ancd\envs\my_pytorch\python.exe` (3.13.15, conda-forge)
- Key packages: numpy, pandas, scipy, scikit-learn, statsmodels, patsy, PyYAML, joblib, matplotlib,
  ViennaRNA 2.7.2 (import name `RNA`). Full audit: `metadata/environment.txt`.
- Recreate the pinned environment with `conda env create -f environment.yml`.
- No GPU required. ~16 GB RAM is sufficient.

## Run order

Each stage is a standalone executable that reads its frozen YAML config and writes deterministic outputs.
Stages must be run in order; stages 09/09a/09b/09c compose the HeLa replication.

| Order | Command | Stage | Outputs (examples) |
|---|---|---:|---|
| 1 | `python src/03_audit_inputs.py` | input audit | `metadata/audit/03_*` |
| 2a | `python src/02_check_reference_compatibility.py data/raw_processed/GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities.txt.gz data/reference/ensembl/release-74_GRCh37/Homo_sapiens.GRCh37.74.gtf.gz` | HEK293T reference check | stdout audit |
| 2b | `python src/02_check_reference_compatibility.py data/raw_processed/GSE145805_RAW.tar data/reference/ensembl/release-88_GRCh38/Homo_sapiens.GRCh38.88.gtf.gz --member GSM4333258_HeLa.out.txt.gz` | HeLa reference check | stdout audit |
| 3 | `python src/04_validate_coordinate_mapping.py` | coordinate check | `metadata/coordinate/04_*` |
| 4 | `python src/05_compute_initial_intersection.py` | HEK293T intersection | `data/interim/05_*` |
| 5 | `python src/06_build_analysis_dataset.py` | HEK293T dataset | `data/final/06_*`, `results/06_*` |
| 6 | `python src/07_run_association_analysis.py` | HEK293T association | `results/tables/07_*` |
| 7 | `python src/07_plot_association_results.py` | figure 2–3 | `results/figures/07_*` |
| 8 | `python src/08a_compute_predicted_structure.py` | predicted structure | `data/final/08a_*` |
| 9 | `python src/08_run_prediction_modeling.py` | HEK293T models | `results/tables/08_*`, `results/models/08_*` |
| 10 | `python src/08b_run_structure_diagnostics.py` | structure diagnostics | `results/tables/08b_*` |
| 11 | `python src/09_build_hela_dataset.py` | HeLa dataset | `data/final/09_*`, `results/09_*` |
| 12 | `python src/09a_compute_hela_predicted_structure.py` | HeLa predicted structure | `data/final/09a_*` |
| 13 | `python src/09b_run_hela_association.py` | HeLa association | `results/tables/09b_*` |
| 14 | `python src/09c_run_hela_model_transfer.py` | HeLa transfer | `results/tables/09c_*` |
| 15 | `python src/10_generate_report_figures.py` | figures 1 & 5 | `results/figures/10_*` |
| 16 | `python src/11_verify_release.py --write-manifest` | freeze release hashes | `metadata/release_manifest.csv` |

Run the test suite from the project root with:

```text
E:\ancd\envs\my_pytorch\python.exe -m pytest tests -q
E:\ancd\envs\my_pytorch\python.exe src/11_verify_release.py --verify
```

## Inputs

Immutable downloads live in `data/raw_processed/`; hashes and provenance are in
`metadata/source_manifest.csv` and `metadata/download_checksums.csv`. The HeLa icSHAPE archive
`GSE145805_RAW.tar` is read in place (member `GSM4333258_HeLa.out.txt.gz`) — no extraction file is
created.

## Key decisions

- Response: count-weighted `combined_ratio` = (Acov1+Acov2)/(AGcov1+AGcov2).
- Main structure metric: `R_flank10` = mean valid reactivity over −10..−1 and +1..+10, each side ≥70%
  valid.
- Grouping / split: connected components of `analysis_gene` ∪ identical 201-nt sequence; 80/20 frozen.
- Primary model comparison: ΔMAE = MAE(M2) − MAE(M4) (HEK293T); external comparison MAE(M1) − MAE(M3)
  (HeLa transfer).
- Frozen deviations and gate decisions: `metadata/08_prediction_modeling_decision_log.md` and
  `metadata/09_hela_replication_decision_log.md`.
- Ridge alpha grids must bracket every selected optimum; stages 08 and 08B fail rather than freeze a
  model when the selected alpha lies on either grid boundary.

## Release integrity

`metadata/release_manifest.csv` records SHA-256 and byte size for code, configs, tests, final/interim
derived data, models, tables, reports, and figures. Raw downloads and references are covered separately
by `metadata/download_checksums.csv` and `metadata/02_reference_manifest.csv`.

## Reports

- Chinese final report: `report/ModStruct_report_zh.md`
- Defense outline: `report/ModStruct_defense_outline_zh.md`
- Per-stage technical reports: `results/*_report.md`

## Note on write permissions

The analysis scripts write under the project workspace. When executed through the DSH harness, subprocess
writes may be denied in the default sandbox mode; run the scripts with the workspace-write / full-access
permission the environment requires (see the sandbox escalation prompt) before the stage can create its
`results/`, `data/final/`, and `data/interim/` outputs.
