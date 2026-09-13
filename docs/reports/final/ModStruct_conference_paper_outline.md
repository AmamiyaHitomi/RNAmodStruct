#ModStruct Conference Paper Outline

## Paper positioning

**Recommended Title**

> **ModStruct: A Reproducible Framework for Quantifying Experimental RNA Structure–m6A Associations and Predictive Value**

**Chinese working title**

> **ModStruct: a reproducible analytical framework for quantifying experimental RNA structure-m6A association and predictive value**

**Paper type**: Empirical research in computational biology, based on a repeatable analysis framework, and does not claim to propose a new prediction algorithm.

**Applicable format**: Full text of a general bioinformatics or computational biology conference; if the target conference has not yet been determined, plan for approximately 4,500–5,500 English words and 6–8 pages of text. Abstract, page, figure, and reference limits should be adjusted after the final meeting.

**Core Research Question**

> How robust is the correlation between local experimental RNA structure and modification levels in the m6A sites that have been called and passed quality control; and can the experimental structure provide transferable prediction gains beyond sequence and background variables?

**One sentence conclusion**

> ModStruct identified small, directionally consistent positive correlations between local experimental structures and m6A levels in HEK293T and HeLa, but no consistent cross-cell line prediction gain was observed for experimental structures beyond sequence and background variables.

## Contribution Statement

The main text should be organized around the following three contributions:

1. **Repeatable data integration framework**: Verify public data sources and field semantics, unify cross-reference versions, cross-transcripts, and link coordinates, and fully record site screening and deletions.
2. **Clear separation of association and predictive value**: Use gene clustering inference to quantify the association of experimental structures with m6A levels, while distinguishing the predictive contributions of sequence, computational structures, and experimental structures through nested ablation.
3. **Strong Cross-Cell Line Validation**: Developed in HEK293T, independent review of associations and direct migration of frozen models in HeLa while checking for gene and identical sequence overlap.

Contents that cannot be included in the contribution statement: first discovery that m6A is related to RNA structure, establishing a universal structural grammar, proving that structure affects m6A, proving that m6A causes structural changes, and proposing a predictive model that is better than SMART-m6A.

## Glossary

| Concept | Uniform text writing | Explanation boundaries |
|---|---|---|
| Projects/Frameworks | ModStruct | Repeatable analysis framework, not called a new deep learning model |
| Experimental structure | experimental RNA structure / icSHAPE reactivity | Structure-dependent reactivity, not equivalent to base pairing probability |
| Computed structure | sequence-derived predicted structure | Structural characterization of ViennaRNA calculated from sequence |
| Modification outcome | m6A level / m6A modification level | GLORI quantitative ratio, range 0–1 |
| Main structure indicator | local experimental reactivity (`R_flank10`) | Mean effective reactivity 10 nt either side of the center |
| Comprehensive effect | pooled association estimate | β/SD=0.01007, about 1.01 percentage points/SD |
| Prediction gain | incremental predictive value | The error difference on the same test population before and after adding the structure |
| External review | cross-cell-line validation / frozen-model transfer | HEK293T is tested directly at HeLa after development |

## Abstract (200–250 words)

### Background and gaps (2–3 sentences)

- There are known links between RNA structure and m6A deposition and regulation, and recent structure-aware models have shown that sequence-predicted structures can improve some m6A classification tasks.
- It is unclear whether experimentally measured intracellular structure correlates stably with quantitative m6A levels and whether it provides predictive gain across cell lines beyond sequence and background variables.

### Method (3–4 sentences)

- Introducing ModStruct: Integrating GLORI single-base m6A quantitative data with icSHAPE experimental structural data.
- HEK293T for discovery and model development, HeLa for correlation review and frozen model migration.
- Association inference is clustered by genes; predictions employ M0–M4 nested feature ablation and block gene and identical sequence leakage.

### Results (4–5 sentences)

- The main effect is 0.01430/SD for HEK293T and 0.00947/SD for HeLa, in the same direction.
- The combined two-cell line fixed effect estimate is 0.01007/SD (95% CI 0.00713–0.01301), equivalent to approximately a 1.01 percentage point increase in m6A proportion.
- Experimental structures do not steadily improve HEK293T internal predictions.
- After the frozen model is migrated to HeLa, the experimental structure worsens the MAE relative to the sequence baseline by about 0.064 percentage points; the strict subset results are similar.

### Conclusion (1–2 sentences)

- Experimental RNA structures have reproducible, but small, observational correlations with m6A levels.
- Correlation does not automatically translate into transferable incremental predictive value, and the benefit of structure-aware models depends on the task, structure representation, and validation conditions.

## 1. Introduction (about 650–800 words)

### 1.1 Bidirectional relationship between m6A and RNA structure

**Purpose**: To establish the biological background while making it clear that observational studies cannot determine the direction of action.

- Briefly describe the relationship between m6A deposition, RNA folding, and RNA-binding protein accessibility.
- Explain three possibilities: structure affects modification, modification changes structure, sequence or cellular environment affects both.
- Citing icSHAPE, m6A-switch, and studies on the mechanisms of structure-determined m6A deposition.

**Transition**: There are already mechanistic examples showing that the structure is worthy of study, but it cannot answer whether the quantitative relationship across the entire transcriptome can be stably reproduced.

### 1.2 Development of structure-aware m6A prediction and its unresolved issues

**Purpose**: To establish an accurate relationship with work such as SMART-m6A and limit the innovation space of this article.

- Explain that most forecasting models take sequences or sequence-derived computational structures as their primary input.
- SMART-m6A Computational structural characterization has been shown to improve predictions in some classification tasks, especially candidate sites whose sequences are indistinguishable.
- Point out that the representational gain of computational structures is not equivalent to experimental structures providing independent cell state information.
- Raising the Gap: The incremental value of experimental constructs for quantifying modification levels, and whether this value can be transferred across cell lines, needs rigorous testing.

**Transition**: Therefore, the focus of this article is not on building more sophisticated predictors, but on building a fair, reviewable evidence framework.

### 1.3 Objectives and design of this study

**Purpose**: Give the problem, framework and verification strategy in one paragraph.

- Introducing the three levels of ModStruct: data alignment, correlation estimation, and nested prediction ablation.
- Design given for HEK293T discovery/development and HeLa independent review.
- Previews two key findings: a small combined correlation of approximately 0.01/SD and no stable gain in prediction across cell lines.

## 2. Materials and Methods (about 1,300–1,600 words)

### 2.1 Study design and analytical scope

**Purpose**: Define estimation objects and analysis boundaries.

- The study objects are m6A sites that have been detected, supported by double repeats, and have qualified structural coverage.
- Q1: Adjusted correlation of local experimental structure with quantitative m6A levels.
- Q2: The incremental predictive value of experimental structures beyond sequence, background and calculated structures.
- Q3: Can correlation and prediction conclusions be reviewed in HeLa.
- Emphasize that this study does not evaluate whether m6A occurs for all adenosines, nor does it make causal inferences.

### 2.2 Public datasets

**Purpose**: Describe data sources, cell lines, technology, and research roles.

- HEK293T: GSE74353 icSHAPE and GSE210563 GLORI.
- HeLa: GSE145805 icSHAPE versus GSE210563 GLORI normoxic control.
- Give replications, reference versions, data types and the fact that they are not measured together from the same sample.
- Data files, checksums and provenance records point to provenance files.

### 2.3 Coordinate harmonization and quality control

**Purpose**: Highlight the most reusable method parts of ModStruct.

- GLORI field semantics, 1-based coordinates and linkage rules.
- HEK293T transcript position to GRCh37.74/hg19, then liftOver to GRCh38.
- hg38 from HeLa with versioned Ensembl transcript mapping.
- Center base, round trip mapping, sequence integrity, isoform and ambiguity mapping checks.
- Report site losses step by step to avoid interpreting unmatched sites as biologically negative.

### 2.4 Outcome and structural features

**Purpose**: Clarify the biological and statistical meaning of each variable.

- The combined proportion of counts from two GLORI replicates was used as the primary outcome.
- `R_flank10`Definition of , center exclusion rule, and minimum 70% valid value requirement for each side.
- Sequence features: one-hot, 3-mer, DRACH isoform, local GC, etc.
- Compute structure: ViennaRNA unpaired probability, free energy, ensemble diversity and structural entropy.
- Experimental and computational structures must be treated as different pieces of information.

### 2.5 Association analysis

**Purpose**: Provide main statistical models that can be reviewed.

- Model form:`m6A level ~ R_flank10 + prespecified covariates`。
- Covariates: sequence motif, transcript region, local GC, coverage, positional variables and available expression.
- Uses robust covariance clustering by gene and 1,000 gene clustering bootstraps.
- Structural effects are expressed as differences in m6A ratio per 1 SD change in reactivity.
- Sensitivity analysis: center exclusion, complete coverage, replicates analyzed separately, equal weight per gene, no isoform ambiguity and conventional DRACH.

### 2.6 Cross-cell-line synthesis

**Purpose**: Convert the correlation results of the two cell lines into the core comprehensive quantity of the paper.

- Report HEK293T and HeLa effects separately, not just the combined value.
- Fixed effects as the main summary and random effects as the sensitivity results.
- Report β/SD, 95% CI, and I²; illustrates the instability of heterogeneity estimates with only two cell lines.

### 2.7 Nested predictive modeling

**Purpose**: Fairly test the incremental value of different feature blocks.

- M0: common background; M1: background + sequence; M2: M1 + calculated structure; M3: M1 + experimental structure; M4: M1 + two types of structures.
- The main comparison is M3 vs M1; M4 vs M2 is used to examine the experimental structure gain after controlling the calculated structure.
- Use Ridge as a transparent baseline and HistGradientBoosting for nonlinear sensitivity analysis.
- Group genes and connected components of the same sequence, and select hyperparameters within the development set.
- The indicator is mainly MAE, and the difference interval is constructed through paired grouping bootstrap.

### 2.8 Frozen-model transfer to HeLa

**Purpose**: Describe where the most rigorous predictive evidence comes from.

- All preprocessing, feature selection and model parameters were determined using only HEK293T development data.
- HeLa tags are not involved in model selection or refitting.
- Simultaneously report all eligible sites and exclude a strict subset of HEK293T development set gene/sequence overlap.
- Clarified frozen imputation rules for missing expression fields in HeLa.

### 2.9 Reproducibility and statistical safeguards

**Purpose**: Convert engineering work into method credibility rather than listing documents.

- Freeze configuration, random seed, software version, feature list and result hash.
- Distinguish between exploratory analysis and confirmatory external review.
- Disclosure of data ethics and data permission instructions.

## 3. Results (approximately 1,400–1,700 words)

### 3.1 ModStruct harmonizes heterogeneous m6A and structure datasets

**Claim**: ModStruct can transform data from different reference versions and transcript systems into a common auditable unit of analysis.

**evidence**:

- HEK293T: 253,087 GLORI merged sites resulting in 4,409 associated sites and 4,096 model sites.
- HeLa: 164,761 merged sites resulting in 25,996 associated sites and 24,960 model sites.
- Report major losses from structural coverage versus transcript mapping, rather than interpreting non-matches as negative.

**Corresponding display**: Figure 1, Table 1.

**Transition**: After identifying the analyzed population, first test whether structure and modification levels co-vary.

### 3.2 Experimental RNA structure is weakly associated with m6A levels in HEK293T

**Claim**: There is a small positive correlation in HEK293T.

**evidence**:

- n=4,409 loci, 1,324 genes.
- β/SD=0.01430, clustering robust 95% CI 0.00595–0.02264.
- Gene clustering bootstrap 95% CI 0.00650–0.02276.
- The seven preset sensitivity analysis directions are all positive.
- Make it clear that the overall R² cannot be interpreted as the proportion of variation explained by the structure alone.

**Corresponding display**: Figure 2A–C, Table 2.

### 3.3 The association replicates in HeLa and remains small in pooled analysis

**Claim**: The direction of association is replicated across cell lines, but the magnitude of the effect is small.

**evidence**:

- HeLa: n=25,996 sites, 4,829 genes, β/SD=0.00947, 95% CI 0.00633–0.01261.
- Fixed effects combined β/SD=0.01007, 95% CI 0.00713–0.01301.
- Translated to a 1 SD increase in structural reactivity, the m6A proportion increases on average by about 1.01 percentage points.
- I²=11.2%, but only two cell lines, no emphasis on heterogeneity values.

**Corresponding display**: Figure 2D, Table 2.

**Transition**: Reproducible statistical correlations do not guarantee that the structure will improve predictions on unseen samples.

### 3.4 Structural features do not improve held-out prediction in HEK293T

**Claim**: Within the scope of current features and models, the structure does not provide stable internal predictive gain.

**evidence**:

- M1 MAE=0.20546, M2=0.21193, the forecast structure worsens the MAE by about 0.6466 percentage points.
- The error change of M3 relative to M1 is close to zero.
- The improvement in M4 relative to M2 is only 0.0052 percentage points, 95% CI −0.0176 to 0.0267 percentage points.
- Neither the extended structure curve nor the 401 nt window changed the conclusion.
- The HEK293T test set was viewed, so this section serves as internal descriptive evidence.

**Corresponding display**: Figure 3A–B.

### 3.5 Experimental structure provides no incremental gain in frozen transfer to HeLa

**Claim**: The most rigorous external tests found no transferable predictive gain from the experimental structure.

**evidence**:

- All HeLa sites: M1 MAE=0.21608, M3=0.21673.
- ΔMAE=M1−M3=−0.00064, that is, the deterioration is about 0.064 percentage points after adding the experimental structure.
- The strict subset worsens by about 0.063 percentage points, with the interval also well below zero.
- HGB significantly improves the overall prediction, but there is still no gain in the experimental structure for M4 relative to M2 with ΔMAE=−0.001687.

**Corresponding display**: Figure 3C–D, Figure 4.

### 3.6 In vivo structure shows a stronger association than in vitro structure

**Claim**: Structural associations may depend on cellular context.

**evidence**:

- 3,705 HEK293T isosite comparisons.
- In vivo effect 0.010862/SD, in vitro effect 0.002386/SD.
- Paired effect difference is 0.008476, 95% CI 0.001517–0.015189.
- Results support clues to condition dependence but do not identify protein binding, molecular crowding, or other specific mechanisms.

**Corresponding display**: Can be incorporated into Figure 4, and moved to supplementary materials when space is insufficient.

## 4. Discussion (approximately 800–1,000 words)

### 4.1 Principal findings

**Purpose**: Summarize the results in one paragraph without repeating all the numbers.

- ModStruct establishes a framework for cross-dataset alignment, correlation inference, and nested prediction validation.
- The association of local experimental structure with m6A levels remained positive across both cell lines, with a combined magnitude of approximately 1 percentage point/SD.
- This association did not translate into stable predictive gains within or across cell lines.

### 4.2 Relationship to structure-aware m6A prediction

**Purpose**: Explain why this study and SMART-m6A can be established simultaneously.

- SMART-m6A mainly studies modified/unmodified classification and drives deep models with sequence-derived structures.
- ModStruct mainly studies the quantitative modification level of detected sites, and uses the experimental structure as the core incremental variable.
- Classification and strength regression, computational structure and experimental structure, cross-validation and frozen external transfer to answer different questions.
- The contribution of this article is to limit the conditions for the predictive value of structures, not to declare that structures are generally useless.

### 4.3 Why association may not translate into predictive gain

**Purpose**: Give multiple compatible explanations and avoid single attribution.

- The average effect is small and it is difficult to significantly reduce site-level errors.
- Structural information may overlap with sequence and background variables.
- The experimental structure and GLORI are not measured together on the same sample, and there is a state and batch mismatch.
- `R_flank10`Possibly compressing position-specific patterns; nonlinear model results reduce the possibility of "just because Ridge is too simple" but do not rule out more appropriate position-wise representations.

### 4.4 Strengths and limitations

**Benefits**: Field and coordinate auditing, inference by gene, leak-proof partitioning, freezing external migration, simultaneous reporting of positive and zero-gain results.

**limit**:

- Only the sites that have been detected and have structural coverage are analyzed, and there is a risk of selection and collider bias.
- The two measurements are from different studies and samples and cannot be interpreted causally.
- HEK293T structure coverage is low and screening may not be random.
- With only two cell lines, broad extrapolation of the combined effects is limited.
- No reliable background of all low-modification candidates was obtained, so transcriptome-wide site classification was not completed.
- The current results cannot rule out that other structural features or more suitable models produce gains in specific subgroups.

### 4.5 Implications and next steps

**Purpose**: Convert constraints into testable next steps without expanding conclusions.

- Predefine sequence model low-confidence sites in development data, and then test experimental structural gains in independent data.
- Prioritize the use of combined structure and modification measurements from the same sample.
- Validate regions, motifs, and in vivo/in vitro interactions in more cellular states.
- Multi-modification structural grammar is retained as a long-term assumption.

## 5. Conclusion (about 120–160 words)

- Reaffirm that ModStruct provides an auditable, reproducible structure-modification correlation and predictive value assessment framework.
- Core values ​​given: Comprehensive correlation approximately 0.01/SD.
- Core bounds given: The experimental structure does not predict gain stably across cell lines in the current task and validation design.
- End by emphasizing that evidence of association, predictive value and causal explanation must be assessed separately.

## Main image design

### Figure 1. ModStruct workflow and analytical populations

- A: GLORI and icSHAPE data sources and discovery/review roles.
- B: Coordinate unification, link verification and structure window extraction.
- C: HEK293T and HeLa site screening funnel.
- D: Analysis routes for association, nested ablation, and frozen migration.
- Existing materials:`results/figures/10_data_overview_figure1.*`。

### Figure 2. Reproducible association between experimental RNA structure and m6A levels

- A: HEK293T modified proportional layered local structural curve.
- B: HEK293T main effect and sensitivity analysis forest plot.
- C: HeLa corresponding results.
- D: Forest plot of two cell lines and combined effects, highlighting 0.01007/SD.
- Existing materials:`07_hek293t_association_overview.*`and`10_hela_replication_transfer_figure5.*`, need to be reassembled into a unified multi-panel diagram.

### Figure 3. Nested ablation separates association from predictive value

- A: Information block representation of M0–M4.
- B: HEK293T each model MAE.
- C: Paired ΔMAE and confidence intervals for M3−M1 and M4−M2.
- D: 401 nt and nonlinear model sensitivity results.
- Existing materials:`08_hek293t_prediction_ablation.*`;The sensitivity panel needs to be added based on the results table.

### Figure 4. Frozen cross-cell-line validation

- A: HEK293T training to HeLa freeze transfer diagram.
- B: MAE of the full set and strict subset of HeLa.
- C: Experimental structural increment ΔMAE.
- D: The in vivo/in vitro effects are poor. If the space is insufficient, it will be moved to supplementary materials.
- Existing materials:`10_hela_replication_transfer_figure5.*`。

## Main table design

### Table 1. Datasets, coordinate systems, and analysis populations

- Data set, cell line, technology, reference version, replicates, number of raw sites, number of associated samples, and number of predicted samples.

### Table 2. Primary association and prediction results

- HEK293T, HeLa and combined effects.
- Internal vs. external primary ΔMAE.
- Keep only the estimators and intervals necessary to support the main conclusion.

## Supplementary materials

- Supplementary Methods: field semantics, coordinate transformations, isomer rules, feature lists and software versions.
- Supplementary Figure S1: Coordinate mapping and center base quality control.
- Supplementary Figure S2: Consistency and coverage relationship between two GLORI repetitions.
- Supplementary Figure S3: Seven correlation sensitivity analysis and regional exploration.
- Supplementary Figure S4: 401 nt window with extended experimental structure curves.
- Supplementary Figure S5: Nonlinear model and feature replacement importance.
- Supplementary Figure S6: In vivo/in vitro structural comparison; will not be repeated if the text is retained.
- Supplementary Table S1: Complete data loss and reasons for non-matching.
- Supplementary Table S2: All models, hyperparameters and grouped cross-validation results.
- Supplementary Table S3: HeLa unseen classes, imputation and strict subset auditing.
- Supplementary Note: The results of m5C, m7G, Nm and STM2457 are only used as exploratory extensions and do not contribute to the main conclusion of this article.

## Evidence mapping

| Thesis Claim | Primary Evidence Documents | Location |
|---|---|---|
| Data and coordinates are auditable |`metadata/provenance/source_verification.md`；`metadata/provenance/01_glori_field_semantics.md` | Methods 2.2–2.3；Result 3.1 |
| HEK293T small positive correlation |`results/reports/07_hek293t_association_validation_report.md` | Result 3.2 |
| HeLa correlation recurrence |`results/reports/09b_hela_association_report.md` | Result 3.3 |
| The comprehensive effect is about 0.01/SD |`results/tables/12_meta_analysis.csv` | Result 3.3；Abstract；Conclusion |
| HEK293T unseen structure prediction gain |`results/reports/08_hek293t_prediction_validation_report.md` | Result 3.4 |
| No gain in HeLa frozen migration |`results/reports/09c_hela_model_transfer_report.md` | Result 3.5 |
| Conclusion of structural gain not changed by nonlinear model |`results/reports/17_nonlinear_models_report.md` | Result 3.5；Discussion 4.3 |
| The association is stronger in vivo than in vitro |`results/reports/19_invivo_invitro_report.md`| Result 3.6 or Supplement |
| Research Boundaries and Statistical Risk |`results/reports/phase1_summary_report.md`；`docs/reports/progress/ModStruct_Stage_Research_Report.md` | Discussion 4.4 |

## Writing order

1. Complete Results 3.1–3.5 first, and lock Figure 1–4 and Table 1–2.
2. Write the Introduction again, so that the questions raised and the results correspond one by one.
3. Write Conclusion, fixing the strongest one to support the conclusion.
4. Write a Discussion to explain the relationship with SMART-m6A, predicted no gain and research boundaries.
5. Finally write Methods and Abstract, and complete reference verification and terminology consistency check.

## Minimum enhancements before submission

The first draft of the meeting can be started directly. One of the most worthwhile additional analyzes before submission is to pre-define low-confidence sites using sequence models in the HEK293T development data, and then freeze the definitions to test the incremental predictive value of the experimental structures in the corresponding subset of HeLa. This analysis responds directly to SMART-m6A's "sequence ambiguous candidate" results, but must be reported in full whether positive or negative.

## Current status

This document represents the structure of the paper pending confirmation by the authors and does not imply that target conferences, formats, or final contribution statements have been determined. Once the conference is identified, space needs to be reallocated according to official page count, abstract structure, anonymity requirements, references and supplementary material rules.
