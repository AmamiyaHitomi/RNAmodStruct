# ModStruct Restart Scope List: Fourth Phase Implementation Version

> This article is directly revised from the original `restart_scope_v1`. The 24 items in the original list all retain traceability numbers, but are no longer classified according to "reason for cancellation".
> Instead, it is reorganized into four phases based on research dependencies. The frozen HEK293T–HeLa mainline is the baseline for subsequent expansions and will not be written back or redefined due to reboots.

## Document information (Material Passport)

- Origin Skill: academic-research-suite
- Origin Mode: experiment-plan / scope-revision
- Origin Date: 2026-09-11
- Revised Date: 2026-09-11
- Verification Status: LOCAL_AUDITED
- Version Label: restart_scope_v2_four_phases
- Supersedes: restart_scope_v1
- Evidence boundary: This revision checks the code, configuration, results and data directory of the current workspace; the latest availability of external data still needs to be reviewed online when the second, third and fourth phases are launched.

## 1. Current baseline and restart discipline

### 1. Baseline completed

The original item 24 "HeLa Q3" has been completed and is no longer listed as to-do:

- HEK293T discovery set: 4,409 association analysis sites, 4,096 model sites;
- HeLa review set: 25,996 association analysis sites, 24,960 model sites;
- HeLa correlation direction is consistent with HEK293T;
- The frozen HEK293T model has completed HeLa full set and strict subset migration;
- The current conclusion is that "there is a small positive correlation between local experimental structure and m6A ratio, but no additional predictive gain is observed";
- Mainline versions, environments, testing and release hashes have been frozen.

### 2. Four execution disciplines

1. **No write-back on the main line**: New work adopts new stage numbers, configurations, tables and reports, and does not cover the frozen conclusions of stages 06-10.
2. **Exploratory Transparency**: Restart analysis occurs after the main results are known, and is marked exploratory by default;Only new independent data verification can upgrade the level of evidence.
3. **Audit first and then analyze**: Any external data must first be audited on source, field, reference version, coordinates, link direction, quantitative semantics and integrity.
4. **Setting thresholds phase by phase**: Enter the next phase after reaching the exit threshold in the previous phase to avoid rolling out 24 interdependent directions at the same time.

---

## Issue 1: Method extension for existing data

### Target

Answer questions about the robustness, conditional dependence, cross-cell line aggregation and model boundaries of existing results without adding new external data or changing the frozen main line.

### Include items in the original list

- #2 Identifiable degradation routes for high-dimensional conditional mutual information;
- #6 Formal meta-analysis;
- #16 401 nt window sensitivity;
- #17 Isomer Perception Analysis;
- #18 HistGradientBoosting;
- #19 SHAP / replacement importance;
- #23 Structural entropy `H_pair` independent analysis;
- #24 is regarded as a completed baseline and will not be executed again.

### Implementation package

#### 1A. Summary of dual cell line effects (formerly #6)

- Use the frozen per-SD correlation coefficients and standard errors of HEK293T and HeLa respectively;
- Report fixed and random effects summary, heterogeneity `Q`, `I²` and `τ²`;
- When there are only two cell lines, the heterogeneity estimate is unstable, and the per-dataset effect must be preserved at the same time, and the pooled estimate is not regarded as the final true value;
- Output suggestions: `config/12_meta_analysis.yaml`, `src/12_run_meta_analysis.py`, `results/tables/12_*`, `results/reports/12_*_report.md`.

#### 1B. Conditional dependency testing (original #2)

- does not directly claim to estimate high-dimensional `I(M;S|sequence)`;
- Generalized Covariance Measure or equivalent residual-on-residual test using cross-fitting;
- Fit `E[M|X]` and `E[R|X]` respectively, where `X` uses only sequence, region, coverage and preset background covariates;
- All fittings and parameter adjustments are grouped by genes and the same sequence connectivity,Infer clustering or grouping permutation by genes;
- HEK293T is a developmental analysis, and HeLa uses a frozen process for external review;
- Output suggestions: `config/13_conditional_dependence.yaml`, `src/13_run_conditional_dependence.py`.

#### 1C. Structural entropy independence assumption (original #23)

- Use the existing `mean_pairing_state_entropy` to examine the association of structural uncertainty with the m6A ratio;
- Simultaneously report unadjusted, covariate-adjusted, and conditional associations after adding `R_flank10`;
- Control multiple comparisons with other predicted structural features and do not interpret a single significant result as a mechanism;
- After defining the analysis in HEK293T, freeze the same model to HeLa.

#### 1D. Isoform Sensitivity (formerly #17)

- Starting from `06_hek293t_all_isoform_mappings.csv.gz` and `09_hela_all_isoform_mappings.csv.gz`;
- Compare at least: four estimates of deterministic representation of isomers, single isomer sites, isomer equal weighting, weighting by available structural coverage;
- Multiple isoforms of the same genomic locus must be kept in the same resampling/grouping unit;
- Short read data are not interpreted as true isoform abundances.

#### 1E. 401 nt window sensitivity (formerly #16)

- Re-extract the full 401 nt window from the fixed transcript sequence, without filling in insufficient windows;
- Simultaneous expansion of sequence one-hot, k-mer and ViennaRNA features cannot just expand the experimental structure while maintaining the 201 nt sequence baseline;
- Following the frozen gene-identical sequence connected grouping principle, the 401 nt sequence hash is reconstructed;
- Reported as a sensitivity analysis without replacing the 201 nt master model.

#### 1F. Nonlinear models and explanations (formerly #18, #19)

- First compare Ridge with `HistGradientBoostingRegressor`, and only use development set group cross-validation for parameter adjustment;
- Preserve the nested ablation structure of M0–M4 and cannot only report the optimal model;
- Explain the importance of prioritizing group permutation and retraining ablation; SHAP is diagnostic only,and is clearly not a causal explanation;
- If the tree model does not have a stable external gain, it will be reported as a negative result and the hyperparameter search will not be expanded to chase the test set performance.

### The order of the first issue

`1A meta-analysis → 1B conditional dependence → 1C structural entropy → 1D isomer → 1E 401 nt → 1F tree model/explanation`

### The first phase exit threshold

- Each implementation package has frozen YAML, executable scripts, machine-readable tables, technical reports and tests;
- All preprocessing is fit within the training fold, and there are no new leaks between HEK293T–HeLa;
-Clearly distinguish between preset mainline, subsequent expansion and external review;
- The results, whether positive or negative, are entered into the summary report.

---

## Issue 2: Newly added public data and cross-technology review

### Target

Increase the independence of structure measurement conditions and m6A measurement technology to determine whether Phase 1 observations are driven by a single structure technology, cell line, or batch.

### Include items in the original list

- #12 HEK293T in-vitro icSHAPE;
- #13 Human RNA MaP (U2OS);
- #14 GSE50676 (PARS);
- #15 m6A-SAC-seq (GSE162356);
- Data auditing section of "Multiple Structured Data Sources" in #4.

### Data priority

1. **#12 in-vitro icSHAPE**: closest to the conditions of existing HEK293T in-vivo and GLORI, with the highest priority;
2. **#15 SAC-seq**: Prioritize evaluation of quantitative fields, background sites, and pairability with HeLa icSHAPE;
3. **#13 Human RNA MaP**: Only enter the correlation analysis after finding the quantitative modification data that matches the cell conditions;
4. **#14 PARS**: Modeled separately as an independent technology and not directly merged with icSHAPE original values.

### Unified audit process

- Record the third-level source and checksum of GEO/paper/document;
- Confirm experimental conditions, cell lines, treatment groups and biological replicates;
- Audit coordinate system, reference version, starting from 0/1 and link direction;
- Determine whether the structure value and modification value are continuous quantitative, semi-quantitative or binary detection;
- Quantify mapping loss, structure loss and selection bias;
-Each data set is estimated independently first, and then summarized across data sets.Direct stacking of primitive values ​​is prohibited.

### Entry threshold for the second phase

- The data are still available from authoritative sources and the license allows for research reuse;
- At least one structural data set can form an interpretable site-level joint table with matching m6A data;
- Does not rely on forcing values from different technologies into the same dimension.

### Second phase exit threshold

- Complete at least one in-vivo/in-vitro control or a cross-modification measurement technology review;
- Give stratified effects for technological heterogeneity rather than just the overall average;
- The new data has completely passed the source and coordinate audit at the same level as the main line.

---

## Issue 3: Multi-modified structural niches and maps

### Target

Expand from a single m6A project to a comparative map of "whether different RNA modifications occupy different structural niches"; the database or website only serves as a delivery form after the map is matured.

### Include items in the original list

- #1 Multiple modification comparison / structural grammar;
- #3 Structural association network;
- #4 RMBase v3.0 and multi-dataset merging;
- #7 lncRNA;
- #10 Database/Website;
- #20 retouch/low retouch background classification;
- #21 log-loss → bits.

### 3A. Multi-modified data reconnaissance

- Candidate modifications: m1A, m5C, m7G, Ψ, ac4C, 2′-O-Me, A-to-I;
- Build evidence cards for each modification: detection technology, resolution, quantitative semantics, cell conditions, reference versions, background collections and licenses;
- Single-base binary annotation and continuous quantitative data analysis; RMBase can only assume the presence/absence layer and cannot pretend to be a modification ratio;
- Comparison of structural niches is initiated only after at least two modifications pass the quality threshold.

### 3B. Comparable feature space

- Unify structural feature definitions but allow different techniques to preserve the dataset hierarchy;
- Preset local reactivity, directionality, predicted unpaired probability, structural entropy, sequence motif and transcript region;
- Hierarchical or hierarchical modeling of cell lines, technologies and coverage to avoid interpreting technology differences as modification differences;
- lncRNA is only added after matching structure and modification measurements are obtained and is not extrapolated from the mRNA results.### 3C. Structural association network (formerly #3)

- The two sides of the bipartite diagram are modification type and structure/sequence/regional characteristics;
- Edge weights are derived from comparable standardized effects and uncertainties, and do not use p-values directly;
- Community discovery is a descriptive tool, and its stability needs to be verified through data set resampling and leave-one-out data set analysis;
- Before multiple modifications are implemented, a "network" with only one modification node, m6A, will not be established.

### 3D. Background classification and information gain (original #20, #21)

- Start Q4 only after obtaining a reliable background of "well measured but low/unmodified";
- Background must match by coverage, detectability, base type, region and expression;
- Classification models report continuous performance beyond calibration, AUROC, AUPRC and decision thresholds;
- `Δbits = [LogLoss(baseline) − LogLoss(augmented)] / ln 2` is only used for qualified probability models;
- If the background still comes from the complement of the list of significant sites, items 20, 21 remain blocked.

### 3E. Database/Website (formerly #10)

- Not a prerequisite for the third phase;
- Only build when at least two modifications, two independent data sets and a stable data dictionary are in place;
- Prioritize publishing downloadable tables, versioned API/schema and static documents, and then build interactive websites last.

### The third phase exit threshold

- Comparable evidence for at least two modifications;
- Modification differences and detection technology differences can be distinguished in the model;
- Each edge in the graph can be traced back to the data set, model and uncertainty;
- Once the website is launched, it must be bound to versioned data publishing instead of just a visual shell.

---

## Issue 4: Disturbance, directionality, causal evidence and conditional advanced modeling

### Target

Advance from observational correlations to directional and intervention evidence; advanced models only serve clear causal or cross-dataset questions, and do not "train large models" as an outcome in themselves.

### Include items in the original list

- #5 Directionality / Causal / Mechanism of the four-layer causal framework;
- #8 Causal Inference (merged with #5);
- #9 Basic model training;
- #11 Deep Learning/Transformer (merged with #9 into conditional methods branch);
- #22 Large-scale hyperparameter search.### 4A. Disturbance data priority

- Prioritize finding writer/eraser perturbations, modification enzyme knockdown/knockout, identical sequence modifications and unmodified controls;
- Both modification changes and structural measurements must be available, and the time sequence, processing conditions, batches and controls must be identifiable;
- Draw causal diagrams prior to analysis to identify treatment, outcome, confounding, mediating and selection variables;
- Do not use ordinary cross-sectional correlations, predictive significance, or SHAP as a proxy for causal effects.

### 4B. Four levels of evidence gradient

1. **Association**: The current main line has been completed;
2. **Directionality**: Time series, before and after disturbance, or asymmetric prediction can only provide directional clues;
3. **Causal evidence**: Requires defensible intervention controls, differential designs, or other identifiable strategies;
4. **Mechanism**: Site-level mutations, structural verification or biochemical experiments are required, and pure calculation results must not be upgraded to mechanistic conclusions.

### 4C. Basic model / Transformer (original #9, #11)

- Not enabled by default; re-evaluated only when a sufficiently large multi-modification and multi-condition corpus is formed in the third period and clear questions that cannot be answered by a simple model are raised;
- Ridge, tree models and existing public model baselines must be set;
- Training, parameter adjustment and external validation are grouped by data set/gene/homologous sequence;
- Evaluation focuses on generalization across technologies, across cell lines, or across modifications rather than a single random partition metric.

### 4D. Hyperparameter search (original #22)

- Only executed within frozen development set or nested cross-validation;
- Search space, budget, stopping rules and main indicators are registered before running;
- Test sets and external review sets must not be used to expand the search space;
- If the performance improvement only comes from post-test parameter adjustment, it must be downgraded to exploratory results.

### Entry threshold for the fourth phase

- There is at least one identifiable perturbation data set, or the third phase corpus is sufficient to raise clear cross-domain basic model questions;
- Ethics, permission and data usage conditions checks have been completed;
- Analysis plan freezes before ending tag is read.

### The fourth phase exit threshold

- Causal conclusions strictly match the identified hypotheses and research design;
- Hierarchical reporting of directionality, causal evidence and mechanism evidence;
- The advanced model must outperform the simple baseline on truly independent external data, otherwise it is retained as a negative result.

---

## Overview of the Four Issues

| Issue | Core issue | Original number | Current data status | Priority ||---|---|---|---|---|
| Issue 1 | Conditional dependence, meta, window, isomer and nonlinear model extension | 2, 6, 16–19, 23; 24 Completed | Data ready | Execute immediately |
| Second period | Review of cross-structural conditions and measurement techniques | 4 (Part), 12–15 | Internet access required for review and download | Parallel reconnaissance possible in the first period |
| Phase 3 | Multi-modification structural niches and maps | 1, 3, 4, 7, 10, 20, 21 | Multiple modifications and large gaps in background data | Start after the second phase |
| Issue 4 | Directionality, intervention, mechanism and conditional advanced modeling | 5, 8, 9, 11, 22 | Dependent perturbation/experiment or large-scale corpus | Conditional priming |

## Original 24-item destination index

| Original number | Destination after revision | Status |
|---:|---|---|
| 1 | Issue 3 3A–3C | More data to be revised |
| 2 | First Issue 1B | Ready for immediate implementation, using GCM downgrade route |
| 3 | The third issue 3C | Wait for at least two modifications to pass the threshold |
| 4 | Second phase data audit + third phase map | Phased implementation |
| 5 | Issue 4A–4B | Merged with #8, data to be perturbed |
| 6 | Issue 1 1A | HEK293T and HeLa already available |
| 7 | Issue 3 3B | lncRNA data to be matched |
| 8 | Issue 4A–4B | Merged with #5 |
| 9 | The fourth issue 4C | Conditional expression, not as the main line in the near future |
| 10 | The third phase of 3E | Start after the map matures |
| 11 | Fourth Issue 4C | Merged with #9 |
| 12 | Second Issue | Priority Audit/Download |
| 13 | Second Issue | Need to match modified data |
| 14 | Issue 2 | Independent technology layering |
| 15 | Issue 2 | Priority candidates for cross-modification measurement technology |
| 16 | First Issue 1E | Implementable, need to re-extract 401 nt |
| 17 | Issue 1 1D | Full isomer mapping table can be used directly |
| 18 | Phase 1 1F | Can be implemented immediately || 19 | Issue 1 1F | Replacement priority takes precedence, SHAP only diagnoses |
| 20 | Issue 3 3D | Blocked by reliable low-retouch background |
| 21 | Issue 3 3D | Eligibility probability label dependent on #20 |
| 22 | The fourth issue of 4D | Freeze development process only |
| 23 | Issue 1 1C | Required features already exist |
| 24 | Baseline completed | No longer on the backlog |

## Next execution entry

The first phase is currently the only phase that does not require new data or experimental authorization and can directly improve the integrity of the paper's evidence. When starting, a new independent
`12–17` phase configuration and script, and complete 1A meta-analysis first; the second phase of data availability reconnaissance can be carried out in parallel, but after the audit
Do not proceed to formal analysis until completed.