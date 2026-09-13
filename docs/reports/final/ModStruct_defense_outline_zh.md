# ModStruct Defense Outline

## Positioning in one sentence

Using the public GLORI (m6A quantification) and icSHAPE (RNA structure detection) data, at the mRNA sites that have been detected and passed quality control,
Answer two-level questions: **whether the local structure signal is related (correlated) to the modification ratio**, and **whether the structural signal can improve the modification ratio
Predictions (predicted value)** and reviewed on HeLa.

## Suggested order of presentation (about 10–12 minutes)

1. **Problems and Motivation (1 min)**
   - m6A is known to be associated with local structures (2015 icSHAPE, m6A-switch; 2024 Shachar; 2026 SMART-m6A).
   - This project does not strive to be the first to "whether it is relevant", but to do a clear and reproducible demonstration: separate "sequence representation gain" and
     "Experimental measurement gain" is the condition for establishing quantitative conclusions under strict division and external review.

2. **Data and Alignment (2 min)**
   - Figure 1: Source, sample conditions, screening process and final number of sites/genes for the two sets of data.
   - Emphasize the difficulty: alignment across studies, years, and reference versions (HEK293T requires GRCh37→hg38 liftOver, HeLa requires
     hg38 direct mapping + versioned Ensembl ID); center base must be A, reverse strand reference is T.

3. **Q1 Association (2 min)**
   - Figure 3 (effect forest plot): R_flank10 corresponds to approximately +1.43 percentage point modification ratio per SD, gene clustering bootstrap
     CI excludes zero; all seven sensitivities are positive.
   - Boundaries: Observability, small effect, 3′UTR dominance, non-linearity.

4. **Q2 Prediction (2 min)**
   - Figure 4 (M0–M4 MAE): The predicted structure M2 is 0.65 percentage points worse than the sequence baseline; the experimental structure M3 is not better than M1
     Reliable increment.
   - Key point: **Statistical correlation ≠ Additional predictive value**; the predictive structure is a deterministic function of the sequence, and the gain can only be called "representation gain".

5. **Q3 Review (2 min)**
   - Figure 5: The correlation is reproduced in HeLa (the direction is consistent, the interval excludes zero, and the amplitude is slightly smaller);However, there is no gain in frozen model migration.
     M0 is the best after migration.
   - Migration failure cannot be solely attributed to "cell specificity" (batch, coverage, culture condition entanglement).

6. **Conclusions and Limitations (1–2 min)**
   - Conclusion: The association between the two cell lines was consistent; no predicted gain was detected.
   - Proactive explanation: Zero gain and migration failure are truthful deliveries, not failed research; Q4 has not been done (lack of background data).

## Expected questions and response points

- **"Is it a contradiction that there is correlation but not prediction?"** No contradiction. Correlation answers "Do structural signals and modification ratios co-vary under given covariates";
  Prediction answers "Can structural signals reduce prediction errors on unseen samples?" The two are different estimators.
- **"Why is the prediction structure useless?"** It is a deterministic nonlinear transformation of the sequence, and the limited Ridge model may not benefit; this is at the representation level
  Conclusion, does not mean that computational structures are biologically irrelevant.
- **"What does HeLa migration failure mean?"** Description The sequence/structural features fitted on HEK293T are not generalized; cell lines, batches, coverage,
  The differences in expression are entangled and cannot be attributed to a single factor.
- **"Is the sample screening loss large and biased?"** Yes. Only a small part of GLORI sites can be mapped to icSHAPE, and it depends on detection and coverage.
  selection/collider bias; the conclusion is limited to the population that has been "detected and passed quality control".
- **"Can we say that structure affects modification?"** No. Observational, cross-study designs cannot determine direction; they can only say "correlated."

## Defense materials list

- Report: `docs/reports/final/ModStruct_report_zh.md`
- Figure 1–5: `results/figures/` (PNG/PDF/SVG)
- Data dictionary and sources: `metadata/provenance/data_dictionary.md`, `metadata/manifests/source/source_manifest.csv`
- Reproduction instructions: root directory `README.md`