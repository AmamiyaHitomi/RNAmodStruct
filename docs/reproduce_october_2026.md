# Reproducing the October 2026 RNAmodStruct analysis

Use this guide with the repository revision associated with the processed-input archive `RNAmodStruct_processed_inputs_oct2026_v1.zip`. The archive SHA-256 is `4da8d0540ed57a84bba5ba64ba220087f67a2dd8bfb1fc03c5fa9452eb308e66`. Keep the ZIP intact. The replay entry point verifies its six input files against the SHA-256 manifest inside the ZIP.

## Inputs

| Input | Purpose |
| --- | --- |
| `06_hek293t_site_level_dataset.csv.gz` and `09_hela_site_level_dataset.csv.gz` | Site-level records used for coverage and selection comparisons |
| `06_hek293t_main_analysis_dataset.csv.gz` and `09_hela_main_analysis_dataset.csv.gz` | Harmonized sites used for the association models and to derive the prediction cohorts |
| `08a_hek293t_predicted_structure_features.csv.gz` and `09a_hela_predicted_structure_features.csv.gz` | Sequence-derived structure features for the historical prediction comparisons |

The [data dictionary](../metadata/provenance/data_dictionary.md), [source verification notes](../metadata/provenance/source_verification.md), [source manifest](../metadata/manifests/source/source_manifest.csv), and [reference manifest](../metadata/manifests/source/02_reference_manifest.csv) describe the input fields, accessions, and reference versions. The primary datasets use GEO GSE210563 (GLORI), GSE74353 (HEK293T icSHAPE), and GSE145805 (HeLa icSHAPE). The HeLa track is the member `GSM4333258_HeLa.out.txt.gz` of `GSE145805_RAW.tar`.

## Run the analyses

Run these commands from the repository root. Change the ZIP path to the downloaded archive. Each replay creates a new output directory, runs the existing analysis scripts there, compares the selected output CSVs with the tracked reference tables, and writes `validation.json` on success. Use a new output path for each rerun.

```powershell
conda env create -f environment.yml
conda activate rnamodstruct

$archive = 'C:\path\to\RNAmodStruct_processed_inputs_oct2026_v1.zip'
Get-FileHash -LiteralPath $archive -Algorithm SHA256
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'

python src/reproduce_october_2026.py replay $archive replay_oct2026_association
python src/reproduce_october_2026.py replay $archive replay_oct2026_prediction --scope prediction
python src/reproduce_october_2026.py replay $archive replay_oct2026_nonlinear --scope nonlinear
python src/reproduce_october_2026.py replay $archive replay_oct2026_matched --scope matched
```

| Scope | Reference tables checked |
| --- | --- |
| `association` | Two primary and two sensitivity tables, plus the two-cell-line synthesis |
| `prediction` | HEK293T test metrics and split audit, plus HeLa transfer metrics and model contrasts |
| `nonlinear` | Historical Stage 17 candidate scores and HeLa prediction metrics |
| `matched` | Supplementary matched-model candidates, selected models, metrics, contrasts, and both site-level prediction tables |

The exploratory paired-condition analysis also requires the two processed reactivity tracks from [GSE74353](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74353). Download them to a directory of your choice and retain the filenames below.

| File | SHA-256 |
| --- | --- |
| `GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities.txt.gz` | `4fd439daad58556b5cc5668feade8d8c0f18e15722927eb2724400c0dfbc47c9` |
| `GSE74353_HS_293T_icSHAPE_InVitro_BaseReactivities.txt.gz` | `8c029a3b1edd904ef8d780580019510f3a274d8466e45c677d0764672a48a58e` |

```powershell
python src/reproduce_october_2026.py replay $archive replay_oct2026_paired --scope paired --geo-directory 'C:\path\to\geo_tracks'
```

The paired scope compares the association and attrition tables. All scopes start from harmonized inputs. To inspect upstream downloads, coordinate mapping, and feature generation, use the manifests above and the numbered workflow in the [project README](../README.md).
