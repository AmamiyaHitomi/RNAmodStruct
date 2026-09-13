# Stage 17: Nonlinear Model and Interpretation

## Material Passport

- Origin: phase-1 exploratory extension
- Verification Status: ANALYZED

Both Ridge and HistGradientBoostingRegressor retain M0–M4 nested ablation, and only the gene-identical sequence connected group cross-validation of HEK293T is used for parameter adjustment; HeLa does not participate in parameter adjustment. M4 external MAE: Ridge=0.21463, HGB=0.18942, HGB’s ΔMAE relative to Ridge=0.025204.

Interpretation of the main results as feature block permutation importance in grouped cross-validation on HEK293T, supplemented by retraining ablation on HeLa. Significance is not a causal explanation; SHAP does not serve as a basis for conclusions, nor does it expand the hyperparameter space based on HeLa performance.
