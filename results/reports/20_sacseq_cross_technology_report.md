# Stage 20: SAC-seq cross-m6A measurement technology review

## Material Passport

- Origin Skill: academic-research-suite
- Origin Mode: experiment-agent/run
- Origin Date: 2026-09-11
- Verification Status: ANALYZED
- Version Label: phase2_stage20_v1

The HeLa master review was performed on 2,440 sites with both valid icSHAPE, complete 201 nt sequence, predicted structure, and SAC-seq quantification. It is a cross-study review of matched cell lines, not the same batch of samples. The 463 sites of HEK293/HEK293T were only subjected to close-match sensitivity analysis.

In the same HeLa site, the adjusted effect of experimental structure on GLORI was 0.004069/SD (95% CI -0.004355–0.012493) and on the SAC-seq calibration score was 0.001210/SD (95% CI -0.004653–0.007073). The Spearman correlation between the two m6A values ​​at these sites was 0.400.

In nested gene-identical sequence group cross-validation, the MAE improvement by adding the experimental structure is M1−M3=-0.000066, and the improvement based on the existing predicted structure is M2−M4=-0.000090. All encoding, scaling, and alpha selections are restricted to the outer training fold.

The 5th and 7th columns of SAC-seq BED are respectively`reported_mean_mutation_ratio_pct`and`reported_calibrated_m6a_fraction_pct`Save, emphasizing that they are the results of the submitter's processing; the main ending is column 7 divided by 100. Method source: https://pmc.ncbi.nlm.nih.gov/articles/PMC9378555/ .
