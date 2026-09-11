# GLORI processed-table field semantics

Checked: 2026-09-10 (Asia/Shanghai)

## Decision summary

- `Sites` is a **1-based genomic coordinate** and identifies the candidate modified adenosine itself.
  On the `+` strand the genomic reference base is `A`; on the `-` strand the genomic reference base is
  `T`, whose transcript-oriented complement is `A`.
- Use `Ratio` as the primary site-level m6A modification proportion in downstream analysis.
- Retain `NormeRatio` only as a prespecified sensitivity variable. It is a deterministic rescaling present
  in the deposited GEO tables, not the modification-level field named by the paper and current README.

## Evidence chain

1. The primary paper defines site methylation level as the A rate: the percentage of reads retaining A in
   total site coverage. It also states that GLORI-tools reports the corresponding A rate as modification
   level: [Liu et al., Nature Biotechnology (2023)](https://doi.org/10.1038/s41587-022-01487-9).
2. GLORI-tools' README calls `Ratio` the A rate or site methylation level and defines `AGcov`, `Acov`, and
   `Genecov`: [GLORI-tools README](https://github.com/liucongcas/GLORI-tools#52-glori-sites-files).
3. The locally pinned GLORI-tools snapshot is commit
   `def608b66ce4a32d18d76b4742df59dd3aaaa5c7` (2026-03-29). The code writes
   `pileupcolumn.pos + 1`, filters transcript-oriented reference A positions, computes the count ratios,
   and performs the statistical test.
4. All 708,103 data rows across the four downloaded GEO tables satisfy the two numerical identities below.
   No row violates either identity at the precision stored in the files.

## Coordinate convention

`pysam` pileup positions are 0-based internally. GLORI-tools writes `pileupcolumn.pos + 1` to the pileup
table (`pileup_genome_multiprocessing.py`, line 42 in the pinned snapshot). The formatter subsequently
passes only `+`-strand positions whose reference base is `A` and `-`-strand positions whose genomic
reference base is `T` (`m6A_pileup_formatter.py`, lines 270-279). The unchanged `pos_1` is propagated to
the final `Sites` column by `m6A_caller.py`.

Therefore a GLORI row is best represented as a one-base, 1-based closed site `(Chr, Sites, Strand)`. If a
BED interval is needed, convert it to `[Sites - 1, Sites)`. A center-base assertion must expect genomic `A`
on `+` and genomic `T` on `-`, while the transcript-oriented center must always be `A`.

## Quantitative fields and units

| Field | Definition | Unit / scale |
|---|---|---|
| `AGcov` | Number of retained reads carrying either A or G at the site after the pipeline's read, base-quality, and A-count filters. | integer read count (effectively UMI-deduplicated molecules in these libraries) |
| `Acov` | Number of those retained reads that still carry A at the site. | integer read count |
| `Ratio` | `round(Acov / AGcov, 5)` in the caller. Every deposited row agrees within the expected five-decimal rounding tolerance. | proportion in `[0,1]`; multiply by 100 for percent |
| `NonCR` | Gene-specific non-conversion background. For transcript-oriented A loci with `AG >= 15`, the formatter accumulates A counts and A+G counts, then computes `sum(A) / sum(A+G)`. Unannotated sites use the configured `ELSE`/overall background. | proportion in `[0,1]` |
| `NormeRatio` | `Ratio * (1 - NonCR)`. This identity is exact for every deposited row. | proportion in `[0,1]` |
| `Genecov` | `sum(A+G counts) / number of qualifying A loci` for the gene; the code calls the denominator `Alength`. It is not total gene reads and not mean depth over every nucleotide. | mean reads per qualifying A locus |

`1 - NonCR` is the corresponding gene-level A-to-G conversion rate. The current public
`m6A_caller_FDRfilter.py` computes a `normedRatio` variable with the formula above but its current output
schema emits `CR = 1 - NonCR` and `Ratio`; the deposited GSE210563 tables instead contain both `NonCR`
and `NormeRatio`. Thus the formula is strongly verified, but the deposited schema reflects the authors'
data-export version rather than the current README table verbatim.

For completeness, the caller's default statistical path performs a one-sided binomial test of `Acov`
successes out of `AGcov` trials against background probability `NonCR`; `P_adjust` is Benjamini-Hochberg
FDR (`fdr_bh`) across the caller output being filtered.

## Downstream field choice

The primary outcome will be `Ratio`, because both the paper and GLORI-tools documentation explicitly
define the site A rate as the m6A level. `NormeRatio` will be retained for sensitivity analysis because it
slightly scales `Ratio` by the estimated conversion rate and is present in the deposited tables. It should
not replace `Ratio` silently, and it should not be described as an independently observed measurement.

Recommended canonical names in derived tables:

- `glori_m6a_ratio = Ratio`
- `glori_m6a_ratio_conversion_scaled = NormeRatio`
- `glori_gene_nonconversion_rate = NonCR`

## Limitations

- The exact historical GLORI-tools commit used to generate GSE210563 was not deposited with each sample.
  The pinned current repository and the deposited data agree on the relevant arithmetic, but the final
  column names differ as noted above.
- `AGcov` and `Acov` are processed read counts, not raw sequencer clusters. Upstream UMI deduplication and
  site/read filters mean they should be treated as effective molecular coverage.
- The gene-level `NonCR` aggregation includes qualifying transcript-oriented A loci as implemented by the
  public formatter; it is a background model input, not a site-specific independent control measurement.

## Local reproducibility pointers

- Code snapshot: `metadata/software/GLORI-tools/`
- Reference assets: `metadata/02_reference_manifest.csv`
- Transcript-version comparison: `metadata/02_reference_compatibility.csv`

AI assistance was used for source retrieval and code/data cross-checking. Every operational conclusion
above is tied to the primary paper, the public pipeline code, or an exhaustive identity check over the
locally deposited tables.
