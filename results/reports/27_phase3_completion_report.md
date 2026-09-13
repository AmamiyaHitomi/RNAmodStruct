# Stage 27: Third Phase Completion and Exit Threshold Report

## Material Passport

- Origin Skill: academic-research-suite
- Origin Mode: experiment-agent/run
- Origin Date: 2026-09-13
- Verification Status: VERIFIED
- Version Label: phase3_stage27_v1

Phases 23–26 of the third phase are all PASS, and the exit thresholds for 6/6 are passed. The official status of Issue 3 is `COMPLETE`; the original multi-modifier dimensions were never merged.

## Main conclusions

- **m5C/HeLa**: n=420, each 1 SD structural change corresponds to a proportional change in m5C -0.003492 (95% CI -0.019900 to 0.012917; FDR q=0.74914).
- **m7G/HeLa**: n=190, interval-level log2 enrichment change -0.026355/structure SD (95% CI -0.187897 to 0.135186; q=0.74914). This result cannot be explained by a single-base m7G effect.
- **Nm/HeLa**: 1281 iso-transcripts, iso-base-paired, poorly structured -0.004505 (cluster bootstrap 95% CI -0.020140 to 0.012171; q=0.74914).
- None of the three main HeLa tests passed FDR 0.05, and the confidence intervals all spanned zero. The safest summary is: **Under the current public data, technology stratification model and coverage threshold, no clear structure-modification association has been detected; this does not prove that the true effect is strictly zero. **

## Keep restrictions

- m7G/HEK293T has only 18 qualified intervals, which are only described and not inferred.
- Nm/HEK is an approximate cell label for HEK versus HEK293T and can only be used as a sensitivity result.
- m5C, m7G, Nm are derived from cross-study observational data; detection efficiency, expression levels and local sequence may still cause residual confounding.
- m1A, Ψ, ac4C and A-to-I stillMaintain HOLD/CONDITIONAL status for technical validity, file size, control reanalysis, or cell matching.

## Next stage

`PHASE3_COMPLETE_GO_TO_PHASE4_PLANNING`: The fourth phase only authorizes entry into perturbation or verification plan planning with controls, which does not mean that there is causal evidence, nor does it automatically authorize large-scale data downloads.