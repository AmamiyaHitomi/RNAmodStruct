# Stage 09 HeLa replication decision and deviation log

## Scope

Stage 09 executes the frozen HeLa replication plan (roadmap W9--W10, gate G5). It uses the same
HEK293T-frozen inclusion rules and the frozen stage-08 models, applied directly to HeLa. No HeLa label,
feature, or parameter was used to select rules or tune models.

## Reference and coordinate decisions

- HeLa icSHAPE is read directly from `GSE145805_RAW.tar` member `GSM4333258_HeLa.out.txt.gz`; the
  extracted member is not materialized in `data/raw_processed/` (immutable download area).
- HeLa icSHAPE transcript IDs are versioned (e.g. `ENST00000606214.5`). They are matched to
  `Homo_sapiens.GRCh38.88.gtf.gz` using `transcript_id "." transcript_version` and to
  `Homo_sapiens.GRCh38.cdna.all.fa.gz` versioned FASTA headers.
- HeLa is on GRCh38, so GLORI sites map directly onto Ensembl release 88 GRCh38 exons; there is no
  hg19→hg38 liftOver step for HeLa (unlike HEK293T).
- 60,350 of 64,616 HeLa icSHAPE transcripts are eligible with exact declared-length agreement across
  icSHAPE, GTF exon sums, and cDNA. Excluded: 3,040 without a GTF transcript and 1,226 without a cDNA
  sequence in release 88.

## Field-semantics deviation (HeLa only)

- The HeLa icSHAPE third field is the literal `*` (no abundance/RPKM is submitted), so
  `icshape_abundance_rpkm` is stored as missing. In the independent HeLa association (09B) the
  `log1p_icshape_abundance_rpkm` covariate is therefore unavailable and excluded, with this deviation
  recorded. In the frozen transfer (09C), the frozen HEK293T encoder imputes that missing covariate with
  the HEK293T development median, which is the correct frozen-preprocessing behavior and is reported as
  a transfer caveat.

## Association replication (09B)

The frozen stage-07 primary model (outcome `combined_ratio`, predictor `reactivity_mean_flank10`,
gene-clustered CR1 covariance, 1,000 cluster-bootstrap resamples) is re-estimated on the HeLa main
dataset. The seven prespecified sensitivity analyses and the nonlinearity diagnostics are repeated
verbatim. This is an independent re-estimate, not evidence that licenses adjusting the HEK293T model.

## Model transfer (09C)

The five frozen HEK293T Ridge models (M0--M4) and their single frozen preprocessing encoder are applied
directly to the HeLa model population. The primary external comparison is ΔMAE = MAE(M1) − MAE(M3).
A strict subset additionally removes HeLa sites whose analysis gene or 201-nt sequence overlaps the
HEK293T development set. Reconstructed frozen feature order is asserted identical to each saved model.

## G5 decision (recorded after 09B/09C)

- **Association replication: PASS.** HeLa re-estimate of the frozen association model is positive and
  excludes zero (+0.00947/SD, cluster-robust CI +0.0063 to +0.0126, p = 3.5e-9), direction-consistent
  with HEK293T (+0.01430/SD), with all seven prespecified sensitivities positive.
- **Predictive transfer: NO reliable increment.** The frozen HEK293T models transfer as-is, but
  experimental structure does not improve on the sequence baseline: ΔMAE = MAE(M1) − MAE(M3) =
  −0.064 percentage points (all) and −0.063 (strict subset), with bootstrap CIs entirely below zero.
  The simple background model M0 transfers with the lowest MAE, indicating the HEK293T-fitted
  sequence/structure features do not generalize their fit to HeLa. Reported as a transfer failure /
  null increment, not attributed to a single cause (coverage, batch, and cell-line differences are
  entangled).
- Strict subset: 22,256 of 24,960 qualifying HeLa sites (2,704 excluded for HEK293T-development
  gene/sequence overlap); results are essentially identical to the all-qualifying results.

These transfer values were regenerated on 2026-09-11 after a development-only alpha-grid audit widened
the HEK293T Ridge grid and selected the interior alpha=300 for M1--M4. No HeLa labels were used for that
correction, and the G5 decision did not change.

These outcomes are recorded in `09b_hela_association_report.md`, `09c_hela_model_transfer_report.md`,
and the corresponding run-status files.
