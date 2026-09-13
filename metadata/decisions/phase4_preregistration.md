# The fourth phase of directional analysis pre-registration and causal diagram contract

## Material Passport

- Origin Skill: academic-research-suite
- Origin Mode: experiment-agent/plan
- Origin Date: 2026-09-12
- Verification Status: LOCAL_AUDITED
- Version Label: phase4_directionality_prereg_v1

## Research questions and main indicators

Is there a pairwise shift in transcript average SHAPE reactivity relative to the DMSO control after 24 h of 40 μM STM2457 treatment in the same HEK293T background? The primary estimator was the median within-transcript difference of STM2457 minus control, using transcript resampling confidence intervals and the Wilcoxon paired test. Gini displacement was the secondary endpoint, and BH correction was performed for both tests.

## Admission and Exclusion

- GSE264642: Same cell line, drug and vehicle control, 2 biological replicates for each probe condition; m6A reduction was verified by meRIP-qPCR in the same experiment, and the access was directional.
- GSE52662 + GSE60034: The structure WT is v6.5 and KO is J1, which violates the same cell background and is excluded from causal statistics.

## Evidence Boundary

STM2457 is not a randomly assigned multi-batch experiment and lacks site-by-site sample m6A quantification, rescue, and inactive compound controls; processing may also affect structure through RNA abundance, RBP occupancy, or off-target pathways. Therefore, the results are at most directional clues and cannot be written as proof of causal effects or mechanisms. The analysis was a retrospective exploratory review of published results.
