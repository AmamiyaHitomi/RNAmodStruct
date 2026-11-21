# Reproducing the supplementary revision analyses

This guide accompanies Supplementary Sections S9.1–S9.3 and the v17
version of Figure S3 from the RNAmodStruct study. Run commands from the
repository root. Use fresh output paths for each rerun.

## Environment and shared inputs

Create and activate the environment in [environment.yml](../environment.yml).
Follow the [October 2026 replay guide](reproduce_october_2026.md) for input
provenance and environment details. Download the
[existing processed-input ZIP](https://github.com/AmamiyaHitomi/RNAmodStruct/releases/download/october-2026-repro-v2/RNAmodStruct_processed_inputs_oct2026_v1.zip).
Its SHA-256 is 4da8d0540ed57a84bba5ba64ba220087f67a2dd8bfb1fc03c5fa9452eb308e66.
Check the digest before extracting it into a clean checkout, preserving the
archive's data/final/ layout.

```powershell
conda env create -f environment.yml
conda activate rnamodstruct
$archive = 'C:\path\to\RNAmodStruct_processed_inputs_oct2026_v1.zip'
$expected = '4da8d0540ed57a84bba5ba64ba220087f67a2dd8bfb1fc03c5fa9452eb308e66'
if ((Get-FileHash -LiteralPath $archive -Algorithm SHA256).Hash -ne $expected) { throw 'ZIP checksum mismatch' }
Expand-Archive -LiteralPath $archive -DestinationPath .
New-Item -ItemType Directory -Path s9_outputs | Out-Null
```

## S9.1 and S9.2

```powershell
python src/p1_coverage_and_within_gene.py s9_outputs/coverage
python src/p1_reactivity_partial_r2.py s9_outputs/partial_R2.csv
python src/p1_local_reactivity.py s9_outputs/positions
python src/plot_p1_coverage_and_splines.py s9_outputs/coverage s9_outputs/coverage_figures
python src/plot_p1_figure_s3_v17.py s9_outputs/positions/positions.csv s9_outputs/figure_s3
```

For figure-only replay, use the archived positions.csv in the command below.

```powershell
python src/plot_p1_figure_s3_v17.py results/tables/revision_s9_v17/positions.csv s9_outputs/figure_s3_from_reference
```

In Figure S3, icSHAPE reactivity is the outcome and the coefficient is scaled
to a 10-percentage-point difference in modification level. Markers use Holm
correction over the 40 position-by-cell-line tests. The plotting script derives
the raw P values for the displayed model from its stored normal-reference
confidence intervals. The p_holm column in the archived positions.csv pertains
to the forward model with modification level as the outcome.

## S9.3

Download the HepG2 and HeLa poly(A) BED files from
[GEO GSE162356](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162356).
Place GSE162356_HeLa_polyA.bed.gz in data/raw_processed/ and pass the
HepG2 BED as a command-line argument. Place GSE145805_RAW.tar from
[GEO GSE145805](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145805)
in data/raw_processed/. The script reads GSM4333262_HepG2.out.txt.gz and
GSM4333258_HeLa.out.txt.gz from that TAR archive. Put the release-88
GRCh38 GTF and cDNA FASTA in data/reference/ensembl/release-88_GRCh38/.
See the [reference manifest](../metadata/manifests/source/02_reference_manifest.csv)
and the [input checksums](../results/tables/revision_s9_v17/sacseq_input_sha256.json)
for exact file names and digests.

```powershell
python src/p1_sacseq_cell_line_check.py data/raw_processed/GSE162356_HepG2_polyA.bed.gz s9_outputs/sacseq
```

## Reference outputs

The [S9 reference tables](../results/tables/revision_s9_v17/) contain summary
results for comparison. The earlier five-scope replay validation applies to
the October 2026 replay guide. These S9 commands describe the additional
revision analyses.