# Stage 15: Isoform Sensitivity

## Material Passport

- Origin: phase-1 exploratory extension
- Verification Status: ANALYZED

All four mapping calibers use genomic sites as analysis/resampling boundaries. Equal weighting and structure coverage weighting are sensitivity weights only; short-read data are not interpreted as true isoform abundances.

- HEK293T / deterministic_representative_isoform：n=4409，β/SD=0.01430（95% CI 0.00595–0.02264）。
- HEK293T / single_isoform_sites：n=2257，β/SD=0.02051（95% CI 0.00881–0.03221）。
- HEK293T / isoform_equal_weight：n=4409，β/SD=0.01623（95% CI 0.00792–0.02455）。
- HEK293T / available_structure_coverage_weight：n=4409，β/SD=0.01623（95% CI 0.00791–0.02454）。
- HeLa / deterministic_representative_isoform：n=25996，β/SD=0.00947（95% CI 0.00633–0.01261）。
- HeLa / single_isoform_sites：n=6050，β/SD=0.01151（95% CI 0.00505–0.01797）。
- HeLa / isoform_equal_weight：n=25996，β/SD=0.00943（95% CI 0.00629–0.01258）。
- HeLa / available_structure_coverage_weight：n=25996，β/SD=0.00944（95% CI 0.00629–0.01258）。
