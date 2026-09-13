# Stage 16: 401 nt window sensitivity

## Material Passport

- Origin: phase-1 exploratory extension
- Verification Status: ANALYZED

The full 401 nt window (unpopulated) retained 2,143 HEK293T and 16,818 HeLa sites. Sequence one-hot, 3-mer, ViennaRNA features reconstructed simultaneously with experimental structure windows; grouping based on genes and identical 401 nt sequence connected components. Hyperparameters are only selected in HEK293T grouped cross-validation, and HeLa only performs frozen external review.

HeLa's M1 MAE=0.22758, M4 MAE=0.22875, ΔMAE(M1−M4)=-0.001164. This result is a 401 nt sensitivity analysis only and does not replace the 201 nt master model.
