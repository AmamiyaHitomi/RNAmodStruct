# Stage 26: Multi-modification hierarchical correlation analysis

## Material Passport

- Origin Skill: academic-research-suite
- Origin Mode: experiment-agent/run
- Origin Date: 2026-09-12
- Verification Status: VERIFIED
- Version Label: phase3_stage26_v1

## Experiment Result

- **ID**: phase3_stage26
- **Type**: analysis
- **Status**: completed
- **Command**: `E:\ancd\envs\my_pytorch\python.exe src\26_run_stratified_multimodification_associations.py`
- **Working Directory**: `E:\my_python\RNAmodStruct`
- **Output Gate**: `GO_TO_PHASE3_SYNTHESIS`

The three HeLa master tests form the same BH-FDR family; HEK Nm is approximate cell line sensitivity, HEK293T m7G is described only because n=18. None of the original modification dimensions were pooled across technologies.

## Hierarchical estimation

- **m5C/HeLa**: mean m5C proportion change -0.003492 (95% CI -0.019900 to 0.012917; q=0.74914; n=420) per ±10 nt mean icSHAPE increase of 1 SD.
- **m7G/HeLa**: log2 peak enrichment change -0.026355 (95% CI -0.187897 to 0.135186; q=0.74914; n=190) per 1 SD increase in intra-peak mean icSHAPE. This is an interval-level association, not a single-base effect.
- **Nm/HeLa**: The structural difference between the modified site and the iso-transcript, iso-base control is -0.004505 (cluster bootstrap 95% CI -0.020140 to 0.012171; q=0.74914; paired n=1281).
- **Nm/HEK**: Approximate cell line sensitivity difference -0.012241 (95% CI -0.027653 to 0.002876; not included in the HeLa master test FDR family).

## Explain boundaries

These are observational associations across studies. The m5C and m7G models adjust coverage/interval width and transcript position respectively, and robustly estimate clustering by gene; Nm control matching reduces transcript and base composition confounding, but cannot exclude residual confounding by expression, detection efficiency, and local sequence environment. Statistical significance does not constitute evidence of causation or mechanism.
