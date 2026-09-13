# ModStruct: Correlation and predictive analysis of human m6A modification levels and local RNA structural signals

**Final report (including HEK293T discovery set and HeLa review set)**

---

## Summary

This study integrates the public GLORI single-base m6A quantitative data and icSHAPE experimental structure detection data.
At the mRNA site, the relationship between the local experimental structure signal (±10 nt average reactivity R_flank10) and the m6A modification ratio was examined, and
Model ablation of sequences, predicted structures, and experimental structures evaluates the additional predictive value of experimental structures, and finally the association and migration are reviewed with HeLa data.

Key results:

1. **Correlation (Q1, HEK293T)**: R_flank10 has a small positive correlation with the m6A ratio. For every 1 standard deviation increase in the structural signal,
   Modification proportion increased by approximately 1.43 percentage points (clustering robust 95% CI 0.60–2.26, gene clustering bootstrap CI 0.65–2.28,
   p = 7.9×10⁻⁴); all seven preset sensitivity analyzes remained positive.
2. **Prediction (Q2, HEK293T)**: On the fixed test set, the ΔMAE of the experimental structure (M3) relative to the sequence baseline (M1) is
   −0.003 percentage points (CI −0.034 to +0.025); the predicted structure (M2) is about 0.65 percentage points worse than M1. not observed
   Reliable additional prediction gain.
3. **Recheck (Q3, HeLa)**: The association was successfully replicated in HeLa (+0.00947/SD, CI +0.0063 to +0.0126, p = 3.5×10⁻⁹,
   The direction is consistent and the amplitude is slightly smaller); however, when the frozen HEK293T prediction model is directly transferred to HeLa, the experimental structure does not bring any gain.
   (ΔMAE = MAE(M1)−MAE(M3) = −0.064 percentage points, strict subset −0.063), and the background model M0 performs best after migration.

Conclusion: **"The local experimental structural signal is positively correlated with the m6A modification ratio" This correlation holds true in both cell lines; however, the experimental structural signal is relatively
Sequential baselines did not exhibit reliable prediction gains. **This conclusion is limited by observational, cross-study, and cross-batch design and cannot be interpreted as causal.

---## 1. Research questions

| Question | Content | Results |
|---|---|---|
| Q1 | Is the modification ratio of the detected m6A sites related to the neighborhood structure signal? | Positive correlation, consistent between the two cell lines |
| Q2 | Do experimental structures improve predictions of modification ratios? | No reliable gain |
| Q3 | Are the results robust across data sets? | Passed the association review; no gain in predicted migration |

The response variable is the continuous m6A modification ratio (the count-weighted combined value `combined_ratio` of GLORI `Ratio`), and the primary structure metric is
`R_flank10` (mean effective reactivity of 10 nt either side of the center, ≥70% effective on either side). The research object is "detected and passed quality control
"site", the conclusion cannot be generalized to the rule for all A sites in the entire transcriptome.

---

## 2. Data and methods

### 2.1 Data

| Purpose | Dataset | Sample | Description |
|---|---|---|---|
| Discovery set structure | GSE74353 | HEK293T icSHAPE in vivo | 10,164 transcripts, GRCh37.74/hg19 |
| Discovery set m6A | GSE210563 | GSM6432590/91 (HEK293T GLORI rep1/2) | GRCh38, processed 14-column table |
| Review set structure | GSE145805 | GSM4333258 HeLa (after NAI-N3 processing) | 64,616 versioned transcripts, hg38 |
| Review set m6A | GSE210563 | GSM6432595/96 (HeLa hypoxia- GLORI rep1/2) | GRCh38, normoxic control |

Coordinate unification: The icSHAPE transcript coordinates of HEK293T were mapped to hg19 through GRCh37.74 annotation, and then aligned through hg19→hg38 liftOver
GLORI's GRCh38 site; HeLa's icSHAPE is hg38 (versioned Ensembl ID), which directly maps to the GRCh38 site and matches
Ensembl release 88. The reverse strand genome reference base is T,But the central base in the direction of the transcript must be A.

### 2.2 Screening process

| Steps | HEK293T | HeLa |
|---|---:|---:|
| GLORI Union | 253,087 | 164,761 |
| Map to icSHAPE transcripts | 20,559 | 86,932 |
| Detected in both replicates | 16,161 | 65,050 |
| Primary structure set (±10 nt ≥70% valid) | 4,409 | 25,996 |
| Model Ready (Full 201 nt) | 4,096 | 24,960 |

### 2.3 Statistics and Models

- **Correlation Model**: `combined_ratio ~ R_flank10 + covariates`, ordinary least squares + robust covariance clustered by genes (CR1),
  and rechecked with 1,000 gene clustering bootstraps. Covariates: DRACH subtype, transcript region, local GC, distance from stop codon,
  Distance from splicing boundary, log1p coverage depth, log1p expression (HEK293T).
- **Prediction model**: M0 (background), M1 (background + sequence one-hot and 3-mer), M2 (+ predicted structure), M3 (+ experimental structure),
  M4 (+both) five-category Ridge regression, alpha ∈ {0.01, 0.1, 1, 10, 30, 100, 300, 1000, 3000}, 50% off GroupKFold
  Adjust parameters; if the optimal value falls on the grid boundary, freeze will be refused. Gene + same 201-nt
  The connected components of the sequence are grouped independently and divided 80%/20% at once.
- **Predicted structure**: Unpaired probability, MFE/nt, ensemble free energy, ViennaRNA 2.7.2 partition function (37°C, dangles=2)
  Pairing state entropy.
- **Review and Migration**: HeLa rebuilds the table using frozen rules and re-evaluates associations independently; the frozen HEK293T model and preprocessing predict HeLa directly,
  Report all eligible sites separately with a strict subset that excludes genes/sequence overlap with the HEK293T development set.

---

## 3. Results

### 3.1 Q1: HEK293T correlation analysis| Indicators | Values |
|---|---|
| locus / gene | 4,409 / 1,324 |
| β (reactivity per SD) | +0.01430 |
| Clustering robust 95% CI | +0.00595 ~ +0.02264 |
| Gene clustering bootstrap 95% CI | +0.00650 ~ +0.02276 |
| p | 7.9×10⁻⁴ |
| Approximately 1.43 percentage points per SD |
| R² | 0.174 |

Seven preset sensitivity analyzes (exclude −2..+2, full coverage only, two GLORI replicates separately, equal weight per gene, no isoform ambiguity only,
Conventional DRACH only) are all forward, BH corrected q ≤ 0.034. Diagnostic shows the presence of heteroscedasticity (inference using clustering robust covariance) and non-linearity
(RESET p = 1.9×10⁻¹⁵, slightly better spline, worse AIC −1.72). In the regional layering, 3′UTR and CDS are forward, and 5′UTR and other regions are
Negative but the interval spans zero.

### 3.2 Q2: HEK293T predicted ablation

| Model | Input | Test MAE (percentage points) |
|---|---:|
| M0 | Background | 20.82 |
| M1 | Background+Sequence | 20.55 |
| M2 | + Predicted Structure | 21.19 |
| M3 | + Experimental Structure | 20.55 |
| M4 | + both | 21.19 |

Main comparison ΔMAE = MAE(M2) − MAE(M4) = +0.005 percentage points (paired groups bootstrap 95% CI −0.018 to +0.027),
No reliable increments of the experimental structure from the "sequence + calculated structure" baseline were detected; the predicted structure (M2) was instead 0.65 percentage points worse than M1. another
There is also no reliable increase in M3 versus M1 of −0.003 percentage points (CI −0.034 to +0.025). 08B development set diagnosis, any extension
None of the experimental building blocks are better than M1, frozen R_flank10; GLORI inter-replicate MAE is 2.65 percentage points (much smaller than the model error of about 20.5),
This shows that repeated noise is not enough to explain the model error.> Note: The HEK293T test set label has been viewed before preprocessing correction, so the internal test results are marked as "reuse test set", descriptive conclusion;
> Independent verification of dependence on HeLa stage.

### 3.3 Q3: HeLa review

**Association review** (25,996 loci / 4,829 genes):

| Indicators | HEK293T | HeLa |
|---|---:|---:|
| β (per SD) | +0.01430 | +0.00947 |
| Clustering robust 95% CI | +0.00595 ~ +0.02264 | +0.00633 ~ +0.01261 |
| bootstrap 95% CI | +0.00650 ~ +0.02276 | +0.00643 ~ +0.01271 |
| p | 7.9×10⁻⁴ | 3.5×10⁻⁹ |

The correlation direction is consistent, the intervals exclude zero, and the amplitude is slightly smaller (about 0.95 vs 1.43 percentage points/SD); all seven sensitivities are positive; the effects are concentrated
At 3′UTR (+0.0124/SD, q = 2.5×10⁻⁹).

**Predicted migration** (frozen model predicts HeLa directly, 24,960 full / 22,256 strict subset):

| Model | Full MAE (percentage points) | Strict subset MAE (percentage points) |
|---|---:|---:|
| M0 | 20.52 | 20.56 |
| M1 | 21.61 | 21.86 |
| M2 | 21.89 | 22.21 |
| M3 | 21.67 | 21.92 |
| M4 | 21.95 | 22.26 |

Main comparison ΔMAE = MAE(M1) − MAE(M3) = −0.064 (all) / −0.063 (strict) percentage points, bootstrap interval is completely negative:
The experimental structure did not bring any gain, and the performance of the frozen sequence/structural features after migration was lower than that of M0 with only background, indicating that the fit on HEK293T
Sequence/structural features are not transferred to their advantage.

---

## 4. Discussion

1. **The correlation is robust and has a small amplitude**. Both cell lines gave positive correlations with intervals excluding zero, with reproducible directions; however, about 1 percentage point per SD
   The magnitude is a "small effect",And there is nonlinearity and regional heterogeneity (dominated by 3′UTR).
2. **Missing predicted gain is a factual result**. Experimental structural reactivity is not a deterministic function of the sequence and should have independent information; but in the current Ridge
   No reliable reduction in prediction error was detected for either HEK293T in-house testing or HeLa migration across baseline and site populations. This is related to
   "There is an association but no predictive gain" is not a contradiction - statistical association and additional predictive value are different estimators.
3. **Migration failure cannot be attributed to a single reason**. HeLa and HEK293T have different research years, processing procedures, and culture conditions. The differences in cell lines and batches
   Differences are entangled; differences in coverage, isoform composition and expression are also candidate explanations and cannot be reduced to "cell specificity".
4. **Explain boundaries**. icSHAPE reactivity is a structural correlation measure, not a pairing probability; observational, cross-study design does not allow for causal inference; adjustment
   Only listed covariates are covered; there are practical boundaries for the treatment of isogenic and homologous sequence dependencies.

## 5. Restrictions

- The research object is "detected and passed quality control" sites, there is selection and collider bias, and it cannot be generalized to the entire transcriptome.
- Only about 52.8% (HeLa) / 8.1% (HEK293T) of the GLORI union sites can be mapped to icSHAPE transcripts, and the screening loss is large and non-random.
- HeLa's icSHAPE does not have an expression field, and the correlation model is missing this covariate; HEK293T is used to develop center digit interpolation during migration.
- The improvement of computational structure is a representation gain and does not mean the discovery of independent biological information outside of the sequence.
- Q4 (modification/low modification classification) has not been carried out due to lack of qualified low modification background data.

## 6. Conclusion

In the two cell lines HEK293T and HeLa, there is a consistent direction and zero interval exclusion between the local experimental structure signal and the m6A modification ratio.
There is a positive correlation; however, the experimental structural signal does not show reliable predictive gain relative to the sequence baseline. This work provides a clear and reproducible
The empirical analysis clearly distinguished the two levels of "correlation" and "predictive value", and truthfully reported zero gain and migration failure.

---

## Appendix

- Figure 1: Data sources and screening process (`results/figures/10_data_overview_figure1.*`)
- Figure 2–3: HEK293T raw trend,Local curve and effect forest plot (`results/figures/07_hek293t_association_overview.*`)
- Figure 4: M0–M4 MAE versus ΔMAE (`results/figures/08_hek293t_prediction_ablation.*`)
- Figure 5: HeLa review and migration (`results/figures/10_hela_replication_transfer_figure5.*`)
- Key tables: `results/tables/07_*`, `08_*`, `09b_*`, `09c_*`
- Instructions for reproducing the operation: see `README.md` in the root directory; see `metadata/` for the source list and data dictionary.