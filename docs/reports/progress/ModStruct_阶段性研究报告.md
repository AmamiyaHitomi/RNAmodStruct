# ModStruct Phased Research Report (Phase 1 to Phase 4)

## Summary

ModStruct studies the relationship between RNA modifications and RNA local structure. The project initially envisaged integrating human transcriptome-wide RNA modification and structure probing data, establishing a map of local structural features of different modifications, and further determining whether these structural features can provide information beyond the primary sequence. Starting from m6A, the research established a complete data audit, coordinate mapping, correlation analysis and predictive modeling process in the HEK293T discovery set and HeLa external review set, and then carried out method expansion, cross-condition and cross-technology review, multi-modification exploration and perturbation directionality analysis in sequence.

By the end of the fourth phase, the project had formed three main insights. First, there is a positive correlation with a more stable direction but smaller amplitude between the m6A ratio and the local experimental structure signal. Second, after information such as sequence, region, coverage, and predicted structure has entered the model, the experimental structure has not yet shown stable external prediction gain. Third, the strength of the association varies with cellular context and measurement technique; exploratory results for m5C, m7G, and Nm are not yet strong enough to support a general “modification-specific structural niche.” Systematic shifts in transcript structural summaries occur after STM2457 treatment, but current evidence only supports directional clues and does not prove that m6A reduction leads to structural changes.

The current four phases of work have completed a research version that can be concluded independently. The most appropriate follow-up line would be to study the conditional dependence of RNA modification-structure associations and their incremental value beyond sequence information, and to validate prespecified local structural patterns using new independent data.

**Keywords:** RNA modification; RNA secondary structure; m6A; icSHAPE; structure detection; interpretable machine learning; conditional dependence

## 1. Research background and issues

The biological functions of RNA are simultaneously affected by primary sequence, spatial folding state, and chemical modifications. Modifications such as m6A, m5C, m7G, and Nm may favor specific sequence and structural environments or may participate in post-transcriptional regulation by altering base pairing, RNA protein binding, or local accessibility. However, published data often come from different cell lines, experimental techniques, reference versions and quantitative systems. If these observations were directly combined, technical differences, coordinate errors, and pseudorepeats could be mistaken for biological regularities.

The project initially proposed the establishment of a "Human RNA modification–structure association atlas",The core hypothesis is that different RNA modifications may occupy different structural niches and form a describable structural grammar. The idea consists of three progressive questions:

1. Are RNA modification levels related to local structural signals?
2. Does the structure still provide independent information after controlling for sequence and other contextual factors?
3. Can observational associations withstand testing across cellular conditions, measurement techniques, and perturbation data?

Research always distinguishes four levels of evidence: Association, Directionality, Causal evidence and Mechanism. Cross-sectional correlations and predictive importance cannot be interpreted as causal effects; differences before and after a disturbance can only provide directional evidence without sufficient comparison.

## 2. Data and overall research design

### 2.1 Baseline study

The project began with m6A to establish a reproducible analytical backbone. HEK293T was used as the development and discovery dataset, and HeLa was used as the external review dataset. Modification data mainly come from GLORI, and experimental structure data mainly come from icSHAPE; predicted structural features are calculated from fixed transcript sequences.

A total of 4,409 HEK293T association analysis sites and 4,096 model sites, and 25,996 HeLa association analysis sites and 24,960 model sites were obtained in the baseline. The primary structural metric is the average experimental reactivity `R_flank10` for 10 nt on either side of the modification site. The data is divided into connected groups of genes and the same sequence as units to avoid the same gene or the same local sequence entering the training set and the test set at the same time.

The analysis is divided into two complementary lines:

- Correlation analysis estimates how much the m6A ratio changes with the structural signal and reports the uncertainty.
- Predictive analysis compares a model using only sequence and background information with a model that adds experimental structural information. Model performance is evaluated using error changes on independent data.

### 2.2 Data quality control

All external data are verified for source, reference genome, coordinate starting method, link direction, transcript mapping, field semantics, file completeness and structural coverage before analysis. Modification values ​​generated by different techniques remain hierarchical, without merging modification ratios, peak enrichments, and binary site states into the same dimension. Multiple nucleotide or transcript observations are also not treated as equivalent independent biological replicates.

### 2.3 The meaning of statistics in this report

To avoid directly equating “a small number”, “a large sample size”, or “a significant p-value” with biological importance,This report interprets statistics as follows:

| Statistics | Meaning in this project | What to pay attention to when reading |
|---|---|---|
| `n` / sites / transcripts | The number of sites or transcripts that actually entered the analysis | Sites and transcripts are the objects being measured and are not equal to the number of independent biological replicates |
| β/SD | How much does the outcome change on average when the structural index increases by 1 standard deviation | When the outcome is a modification ratio of 0–1, 0.010 is approximately equal to 1.0 percentage point; it is the adjusted average association and does not mean that a single site must change by 1% |
| 95% CI | The range of effects that is compatible with the data under the current model and sampling assumptions | The interval spans 0 means that the data allows for positive, no effect, and small negative effects at the same time, and the direction cannot be confirmed |
| p-value | Under a zero-effect model, the probability of observing the current or more extreme outcome | The p-value does not represent the "probability that the conclusion is true" nor does it reflect whether the effect is large enough |
| BH q value / FDR | Corrected significance index after testing multiple hypotheses simultaneously | q<0.05 controls the expected false discovery proportion in a set of results judged to be positive, not the error probability of a single conclusion |
| MAE | The average absolute distance between model predictions and true values | When the modification ratio is expressed as 0–1, MAE=0.20 is approximately equal to an average difference of 20 percentage points; the lower the MAE, the better |
| ΔMAE | Baseline model MAE minus enhanced model MAE | A positive number indicates that the error decreases after adding new features, and a negative number indicates that the error increases instead |
| Spearman ρ | The degree of agreement between the two sets of numerical rankings, ranging from -1 to 1 | ρ=0.40 indicates a moderate degree of same-direction ranking relationship, but does not mean that the two technologies give the same absolute value |
| I² | The proportion of effect differences between different data sets in the summary analysis | There are only two cell lines in this project, and the estimate of I² is very unstable and cannot be used alone to prove consistency |

### 2.4 Overview of key data

| Analysis problem | Actual data | How to interpret the values | Supportable conclusions |
|---|---:|---|---|
| Baseline data scale | HEK293T: 4,409 associated sites, 4,096 model sites; HeLa: 25,996 associated sites, 24,960 model sites | "Associated sites" meet the conditions for association analysis;"Model site" also requires that the sequence and each model feature are complete | Both cell lines have large-scale site-level analysis data, but the independent experimental units are far less than the number of sites |
| Dual cell line m6A–structural summary | Fixed effects β/SD=0.01007, 95% CI 0.00713–0.01301 | A 1 SD increase in structural reactivity is associated with an average increase in m6A ratio of approximately 1.01 percentage points; the compatibility range is approximately 0.71–1.30 percentage points | Supporting small observational associations in the same direction |
| Conditional dependence | HEK293T residual correlation 0.0222, p=0.169; HeLa 0.0141, p=0.0264 | After subtracting preset sequence and background information, the remaining co-variation is only of the order of magnitude 0.01–0.02 | The remaining association is weak and the statistical evidence for the two cell lines is inconsistent |
| 401 nt external prediction | HeLa M1 MAE=0.22758, M4 MAE=0.22875 | After adding predicted and experimental structures to M4, the average error is about 0.00117, or about 0.12 percentage points higher than the sequence model | The expanded window does not show structural gain |
| Nonlinear external prediction | HeLa M4: Ridge MAE=0.21463, HGB MAE=0.18942; ΔMAE of HGB M4 relative to M2=-0.001687 | The overall average error of the nonlinear model is about 2.52 percentage points lower than Ridge; but after adding the experimental structure, the HGB error increases by about 0.17 percentage points | More complex models have overall advantages, but the advantages are not brought by the experimental structure. |
| In vivo vs. in vitro structure | n=3,705 for the same site; β/SD=0.01086 in vivo, 0.00239 in vitro; effect difference=0.00848, 95% CI 0.00152–0.01519 | At these sites, the structural association is approximately 0.85 percentage points higher in vivo than in vitro | Supports condition-dependent cues, does not identify specific intracellular causes |
| Cross-m6A technology review | HeLa n=2,440; GLORI β/SD=0.00407; SAC-seq β/SD=0.00121; ρ=0.400 | There is only moderate agreement between the two technologies for site ranking, and the interval of both structural effects spans 0 | The current matching subset does not provide a clear cross-technology structural association review || Multiple modification analysis | m5C n=420; m7G n=190; Nm paired n=1,281; three-term q=0.749 | The estimates from the three analyzes are different and cannot be added directly; the current uncertainty ranges all include 0 | No clear modification-specific structural association was detected, and the true effect cannot be proven to be zero |
| STM2457 perturbation | n=24,898 transcripts; control median 0.23241, treatment group 0.25278; median pairwise shift 0.02384 | 0.02384 is the normalized SHAPE reactivity unit, not the 2.38% m6A change; it describes the movement of the summary structural signal before and after treatment of the eponymous transcript | Supports structural directional changes after treatment and does not prove that m6A reduction is the sole cause |

## 3. Basic research before the first issue

The first phase is not the starting point for the project. Prior to this, the project had completed a full round of research from data semantics, coordinate mapping and sample construction to HEK293T discovery, predictive modeling and HeLa external review. The first to fourth phases are based on this baseline to continue to test robustness, condition dependence, multi-modification differences and perturbation directionality.

### 3.1 Data semantics and reference system verification

The project first confirmed that GLORI's `Sites` are the 1-based genomic coordinates of modified adenosine; the positive strand records genome A, and the negative strand records genome T, but still corresponds to A along the transcription direction. The primary outcome `Ratio=Acov/AGcov`, range 0–1, represents the site m6A ratio; `NormeRatio=Ratio×(1-NonCR)` is only used for sensitivity analysis. A total of 708,103 rows in the four GLORI tables all satisfy the two numerical identities, and the failure rate is 0. This ensures that the "percentage points" mentioned below do correspond to the modification ratio, not the coverage or peak intensity.

Six data files and nine reference files all passed file size, SHA, and full streaming read checks; all six data sets had malformed records and icSHAPE vector length failures of 0, and the candidate processing data was approximately 458.7 MB. HEK293T icSHAPE uses GRCh37.74/hg19, while GLORI uses GRCh38, thus establishing a chain-wise mapping of "transcript position → legacy exon → hg19 → liftOver → GRCh38";HeLa is independently mapped according to its GRCh38 and Ensembl release 88 compatible systems.

### 3.2 Coordinate verification and data filtering

HEK293T coordinate transformation was tested with 30 hierarchical positions covering the positive and negative strands, exon interiors, and splicing boundaries. The 30 central bases are all A, and all the GRCh38 coordinates of the correct chain orientation are obtained and passed the liftOver round-trip check, with a failure of 0; 12 of them match the GLORI site. This shows that the conversion rules are self-consistent within the test range, but does not mean that all modification sites are covered by structural experiments.

The two GLORI repeats contained 214,715 and 215,021 unique sites, respectively, for a total of 176,649 and a combined total of 253,087. The screening process is as follows:

| Steps | Number of sites | Data meaning |
|---|---:|---|
| GLORI merged set | 253,087 | Appears in at least one modification measurement |
| Mapped to icSHAPE | 20,559 | Two experiments can point to the same RNA location |
| Detected in both GLORIs | 16,161 | Combined repeat counts estimated modification ratio |
| Meets ±10 nt structural window | 4,409 | Can enter correlation analysis |
| Full 201 nt model input | 4,096 | Ready for predictive analysis |

230,866 sites were unmatched due to no qualified icSHAPE transcript exon coverage, accounting for 91.22% of the combined set; there were 939 liftOver failures and 722 ambiguous mappings. This reduction mainly reflects the difference in coverage between the two experiments and cannot be interpreted as “no structure” or “no m6A” at the excluded sites.

### 3.3 HEK293T Discovery Analysis

The main analysis rules are pre-fixed: the site must be detected in both replicates; the outcome is `(Acov1+Acov2)/(AGcov1+AGcov2)`; the structural index takes 10 nt on each side of the center and the average reactivity excluding the center, and at least 7/10 bases on each side are valid. The 4,409 associated sites were derived from 1,324 genes, with a maximum of 38 sites per gene, thus inferring clustering by gene. 7,266 of the 20,559 mappable sites corresponded to multiple compatible isomers, up to 19, which was the basis for subsequent isoform sensitivity analysis.The two GLORI replicates had a Pearson correlation of 0.9879, a Spearman correlation of 0.9728, and a median absolute difference of 0.01549, i.e., the two modification proportions typically differed by about 1.55 percentage points. It illustrates that modification measurements are highly reproducible but does not by itself demonstrate structural association.

The main adjusted model resulted in β/SD=0.01430, gene clustering robust 95% CI 0.00595–0.02264, p=0.000787; 1,000 gene clustering bootstrap intervals 0.00650–0.02276. That is, a 1 SD increase in structural reactivity is associated with an average increase in m6A proportion of approximately 1.43 percentage points. The seven prespecified sensitivity analyzes were all positive, ranging from 0.01020–0.02051/SD, q≤0.0343. The overall R²=0.1742 means that all model variables jointly explain about 17.4% of the site differences, and it cannot be said that the structure alone explains 17.4%. Post hoc partitioning suggested that the 3′UTR and CDS were positive, and other regions were estimated imprecisely and were only used as exploratory clues at the time.

### 3.4 Predictive modeling and internal boundaries

The predictive analysis included 4,096 sites and 4,092 unique 201 nt sequences, divided into 3,271 development sites and 825 retained test sites by connected grouping of genes and identical sequences; overlap between genes, identical sequences, and joint groupings was 0 for both groups. M1 uses background and sequence, M2 adds predicted structures, M3 adds experimental structures, and M4 adds both types of structures.

In the test set, M1 MAE=0.20546 and M2=0.21193. After adding the prediction structure, the error increases by about 0.65 percentage points. M4 MAE=0.21188, which is only an improvement of 0.0052 percentage points relative to M2, and the 95% CI is -0.0176 to 0.0267 percentage points. It cannot be confirmed that the experimental structure brings reliable gains; the extended structure curve in the development set is not better than M1. The HEK293T test set was later used for pretreatment corrections and supplementary diagnostics, so these results were development phase evidence and the project subsequently moved to a fully independent HeLa external review.

### 3.5 HeLa External Review

The two GLORI repeats of HeLa have 137,662 and 140,705 unique sites respectively, and the combined number is 164,761; 86,932 structurally mappable sites are retained in sequence.With 65,050 double repeat sites, 25,996 association sites, and 24,960 model sites, the association set covers 4,829 genes.

HeLa main association was β/SD=0.00947, cluster robust 95% CI 0.00633–0.01261, bootstrap interval 0.00643–0.01271, p=3.53×10⁻9. That is, a 1 SD increase in structure is associated with an average increase in the m6A ratio of about 0.95 percentage points; the magnitude is smaller than that of HEK293T, which is about 1.43 percentage points, but in the same direction. The seven preset sensitivity analyzes were all positive (0.0058–0.0115/SD, q ≤ 0.00048), and the post hoc partitioning suggested that the signal was mainly concentrated in the 3′UTR. This constitutes an external recurrence of the direction of association, which still cannot be upgraded to evidence of causation.

The frozen HEK293T model was directly migrated to 24,960 HeLa sites; after excluding all sites that were isogenic or homogeneous to the HEK293T development set, the strict subset had 22,256 sites. In all HeLa, M1 MAE=0.21608 and M3=0.21673. After adding the experimental structure, the error worsens by about 0.064 percentage points; the strict subset also deteriorates by about 0.063 percentage points. HeLa outcomes and parameters were not involved in model selection or refitting, so this is a more rigorous external test.

### 3.6 Early baseline conclusion

Before the first phase, the project had completed the closed link of “clarifying data meaning—verification of cross-version coordinates—quantitative screening—HEK293T discovery—leakage prevention modeling—HeLa external review” and saved frozen configurations, result hashes, and automated checks. The HEK293T correlation analysis failed 0 out of 14 tests and the overview chart passed all 21 quality checks. This proves that the computational process is repeatable and does not mean that biological causation has been proven.

The early conclusion was that there was a reproducible, small, observational positive correlation between local experimental structure and m6A ratio with an amplitude of approximately 1 percentage point per SD; however, after sequence and background information entered the model, the experimental structure did not show stable additional predictive value. The subsequent four issues focus on these two conclusions, examining their dependence on methods, structural conditions, measurement techniques, modification types, and perturbation backgrounds.

## 4. Issue 1: Method extension of existing data

The first phase examines whether baseline conclusions rely on a specific statistical model, window length, isomer selection, or structural feature definition.Two-cell line pooled analysis yielded fixed effect β/SD=0.01007 (95% CI 0.00713–0.01301), DerSimonian–Laird random effect β/SD=0.01027 (95% CI 0.00675–0.01379), I²=11.2%. With only two cell lines, the stability of heterogeneity estimates is limited, and dataset-by-dataset effects remain the focus of interpretation.

Since the m6A outcome is on a 0–1 scale, β/SD=0.01007 translates to approximately 1.01 percentage points. It means that sites that differ by 1 SD in structural reactivity will, on average, differ by about 1 percentage point in their m6A proportions over the average cases compared by the models; it does not mean that the same change will occur at every site, nor does it mean that changing the structure will cause this change.

Conditional dependence analysis uses the residual correlation test of cross-fitting. HEK293T residual correlation was 0.0222, group sign flip test p=0.169; HeLa residual correlation was 0.0141, p=0.02639. This analysis is not a direct estimate of high-dimensional conditional mutual information, nor is the statistical evidence entirely consistent between the two data sets.

Both residual correlations are very close to zero. A HeLa p-value below 0.05 mainly indicates that a weak residual association can be discerned in a larger sample; it does not interpret 0.0141 as having a strong predictive or biological effect. HEK293T did not pass the same threshold, also limiting cross-dataset consistency.

In structural entropy analysis, after adjusting for covariates and `R_flank10`, β/SD=-0.00328 (BH q=0.6873) for HEK293T and β/SD=-0.00460 (BH q=0.01167) for HeLa. Stable independent effects were not supported in the development set, so the small negative result for HeLa was retained as an exploratory idiosyncratic signal and not upgraded to a mechanistic conclusion.

The isomer sensitivity analysis compared the calibers of deterministic representative isomers, single isomer sites, isomer equal weighting, and weighting by structural coverage. The β/SD range was 0.01430–0.02051 for HEK293T and 0.00943–0.01151 for HeLa, with consistent direction of association. This suggests that the results are not entirely driven by a representative transcript selection rule, but that coverage weights cannot be interpreted as true isoform abundances.

The full 401 nt window retained 2,143 HEK293T sites and 16,818 HeLa sites.The MAE in HeLa is 0.22758 for M1 and 0.22875 for M4, and widening the window and adding structural features still does not produce external gain. In nonlinear modeling, the overall error of HistGradientBoosting is lower than that of Ridge, with the HeLa M4 MAE being 0.18942 and 0.21463 respectively; however, the ΔMAE of HistGradientBoosting's M4 relative to M2 is -0.001687, indicating that the experimental structure has not yet stably improved the existing model.

M1 uses background and sequence, M2 adds predicted structures on this basis, M3 adds experimental structures, and M4 contains both predicted and experimental structures. The lower the MAE, the better. Since the modification ratio is expressed as 0–1, HGB’s M4 MAE=0.18942 can be understood as the predicted value deviating from the measured ratio by about 18.9 percentage points on average. M4 has more experimental structures than M2, but ΔMAE is negative, indicating that the new structures do not reduce the error.

The first issue thus strengthens the original conclusion: the small correlation maintains direction under a variety of analytical calibers, but the independent information and predictive gain of the structure are not consistently supported.

## 5. Issue 2: Cross-structural conditions and cross-measurement technology review

The second phase examines whether associations are driven by single structural conditions or by modified measurement techniques.

A total of 3,705 sites were retained in the isosite comparison of HEK293T in vivo and in vitro icSHAPE. The structural effect was 0.010862/SD in vivo (95% CI 0.001990–0.019735) and 0.002386/SD in vitro (95% CI -0.006233–0.011006). The paired gene bootstrap difference in vivo minus in vitro was 0.008476 (95% CI 0.001517–0.015189). This result suggests that the modification-structure association may depend on the intracellular environment, but it cannot be determined whether it is caused by protein binding, molecular crowding or other factors.

A total of 2,440 sites were included in the HeLa cross-m6A measurement technology review. The adjusted effect of experimental structure on GLORI modification ratio was 0.004069/SD (95% CI -0.004355–0.012493) and on SAC-seq calibration score was 0.001210/SD (95% CI -0.004653–0.007073). The Spearman correlation of the two modification measures was 0.400.Model error for SAC-seq outcomes did not improve after adding experimental structure.

Human RNA MaP and GSE50676 PARS remain HOLD due to lack of quantitative m6A data for matching conditions. The overall conclusion of the second phase is that the observed effects are small, in vivo/in vitro direct comparisons provide clues to differences in structural conditions, and the replication support of cross-technology matching analysis is insufficient; it cannot be determined that there are indeed differences in effects between technologies, nor can the existing results be regarded as constant rules across platforms.

## 6. The third issue: Exploration of multi-modified structural features

The third phase extended the analysis to m5C, m7G, and Nm to examine whether different modifications present different structural niches. Since the three types of modifications are expressed using proportion, interval enrichment, and site status respectively, the study is modeled separately by technology and modification, and the overall effect that lacks common dimensions is not calculated.

Including 420 sites in the HeLa m5C analysis, the m5C proportion changed by -0.003492 (95% CI -0.019900 to 0.012917; FDR q=0.74914) per 1 SD increase in local structural signal. The HeLa m7G analysis included 190 bins, with a log2 enrichment change of -0.026355 (95% CI -0.187897 to 0.135186; q=0.74914) for each SD increase in intrapeak structure. This result is an interval-level estimate and cannot be interpreted as a single-base effect. HeLa Nm analysis resulted in 1,281 iso-transcripts, iso-base pairs, and a structural difference of -0.004505 between the modified site and the control (cluster bootstrap 95% CI -0.020140 to 0.012171; q=0.74914).

All three HeLa main tests failed FDR 0.05, with confidence intervals spanning zero. HEK293T m7G has only 18 qualified intervals, and only descriptive results are retained; HEK Nm is an approximate cell label and only sensitivity analysis is performed. m1A, Ψ, ac4C, and A-to-I remain HOLD or CONDITIONAL due to technical validity, file size, control reanalysis, or cell matching issues.

The third phase demonstrated that multi-modification data can enter the same audit and feature extraction framework, but has not yet established evidence that different modifications have repeatable structural fingerprints. The current data is also insufficient to construct a structural association network in which each edge can be stably traced and reviewed.## 7. Issue 4: Disturbance and Directional Evidence

The fourth phase uses STM2457 and vehicle structural data in the same HEK293T background in GSE264642 to analyze whether there are systematic changes in the RNA structure after METTL3 inhibition occurs. The data contains 2 biological replicates for each of the two treatment conditions, but the transcript statistics table summarized by the authors was used during the current period of the project.

A total of 24,898 paired transcripts were analyzed. The mean SHAPE reactivity median shift for STM2457 minus vehicle was 0.0238433, with a transcript bootstrap 95% CI of 0.0230258–0.0244911. The median shift in reactive Gini was -0.00179552 (95% CI -0.00210929 to -0.00154519).

The median average transcript reactivity was 0.23241 for the control group and 0.25278 for the STM2457 group. The median pairwise shift of 0.02384 is calculated by first calculating "treatment minus control" for each homonymous transcript and then taking the median of these differences, and is therefore not equivalent to a simple subtraction of the medians of the two groups. This value uses normalized SHAPE reactivity units and cannot be written as a 2.38% increase or decrease in the m6A ratio. The small decrease in Gini indicates a slightly more even distribution of reactivity within transcripts, but this interpretation is still limited to aggregated metrics.

This result suggests that there is a systematic shift in transcript structure summaries following METTL3 inhibitory treatment, but it cannot be interpreted as an identified causal effect of m6A deletion. The current data lack site-by-site m6A changes for the same sample, rescue, inactive compound controls, and independent review. Drugs may also affect structure through RNA abundance, RNA-binding proteins, or off-target pathways. Nor are large transcript observations a substitute for independent experimental replication.

Another candidate WT/KO data was excluded from causal statistics due to inconsistent cellular background. Base models, Transformers, and large-scale hyperparameter searches were closed by preset thresholds due to lack of multi-modification, multi-condition, and independent external test sets. The fourth period finally reaches the Directionality layer and has not yet reached the Causal evidence or Mechanism layer.

## 8. Comprehensive discussion

The Phase IV study did not support a simple, general, and strong "RNA modification determines local structure" model.Existing evidence shows a small correlation between m6A and local experimental structures, in vivo/in vitro comparisons provide clues to conditional differences, and cross-technical replication support is still insufficient. Adding structure does not stably improve predictions, which may be related to the overlap of its information with sequence and background variables, and may also be affected by measurement errors, feature expression and model limitations; currently, high overlap of information cannot be regarded as a proven explanation.

The original structural grammar of the project can still serve as a long-term research framework, but currently it can only serve as a hypothesis to be tested. The estimation objects and techniques in multi-modification analysis vary greatly, and it is difficult to completely separate modification types and measurement techniques. At this stage, if a modification structure network or database is directly established, it is easy to show differences in data sources as differences in modification categories.

The value of this study is mainly reflected in three aspects. First, a repeatable process is established from source verification, coordinate unification, structural feature extraction to group modeling and external review. Second, both positive and negative results are retained, with a clear distinction between statistical association, predictive gain, and causal evidence. Third, it reveals the importance of conditions and measurement techniques in modification-structure research, providing a clear design basis for subsequent confirmatory research.

### 8.1 Overall project idea and logical evaluation

The project has now formed a relatively solid data integration, correlation analysis and prediction verification process, but the scientific main line still needs to be further tightened. The most mature evidence comes from the m6A baseline and its robustness and external review; the multi-modification and perturbation analysis is an extended study and has not yet formed a complete chain of mechanism evidence with the baseline. This section is a comprehensive evaluation of the reported results and does not represent adding new experiments or rerunning all analyses.

**Advantages of analytical logic. ** The preliminary research was carried out in the order of data field verification, reference system and coordinate unification, HEK293T discovery, leakage prevention modeling and HeLa external review. First, ensure that the data can be interpreted, and then examine the scientific issues. The project clearly differentiates between "whether structure and modification ratios are relevant" and "whether structure improves predictions outside of sequence and context." Small correlations are not inconsistent with the absence of stable prediction gains; windows, isomers, nonlinear models, and cross-condition checks further define the scope of the results. It is an important advantage of the existing logic that studies retain results that are inconsistent with original expectations.

**Original question versus actual estimated objects. ** The original idea focused on whether different modifications prefer different structural environments; the current most well-documented m6A analysis focuses on the relationship between local structure and modification ratio at sites where modifications have been detected and have qualified structural coverage. The selection of modification sites is a different issue from the modification level of the detected sites, and the latter cannot replace the former. The existing data are more suitable to support the study of correlation between modification level and local structure;To interpret whether any candidate site is modified, reliable unmodified or low-modification background and corresponding design are also required. Therefore, structural grammar remains a hypothesis to be tested.

**Evidence roles at each stage. ** The early baseline and first period constitute the body of correlation, robustness, and predictive gain; the second period examines cross-condition and cross-technology recurrence; the third period evaluates the feasibility of multi-modification expansion; and the fourth period provides clues to structural changes after perturbations. These stages are not progressively stronger proofs of the same cause and effect. Modification proportions, interval enrichment, and matching structure differences in phase III cannot be directly used to rank the structural preferences of modification classes. The fourth phase turned to "whether the structure changes after inhibition treatment", but the lack of modification changes of the same sample and sufficient repeated layer verification cannot yet explain the previous observational correlation. The completion of the fourth phase means that the implementation scope of this round is completed, but does not mean that the closed loop of the mechanism is completed.

**Sample selection and extrapolation boundaries. ** HEK293T was screened from 253,087 merged sites to 4,409 associated sites, and conclusions were first applied to analyzed samples that met modification detection, structural coverage, and other quality requirements. This screen has a basis for quality control, but the differences between selected and unselected sites still need to be evaluated; the results cannot be automatically generalized to all RNAs or all candidate modification sites.

**Conditional differences and cross-technology recurrence need to be explained separately. ** In vivo/in vitro assays directly test for differences in effects and provide clues to condition dependence. GLORI/SAC-seq matching analysis did not confirm a clear association. First of all, it indicates that there is insufficient support for cross-technical replication, and it cannot be concluded that the effects between technologies are different based on the significance status of different analyses. Statistically, "one is significant and the other is not significant" does not mean that the difference between the two effects is significant.

**Core issues that suggest convergence. ** How robust is the association of local RNA structure with m6A levels? Can it provide transferable prediction gains beyond existing sequence and background information? Current conclusions should be qualified as follows: no stable additional predictive gain is seen within the range of data, structural features and models evaluated. It cannot be concluded from this that the structure has no biological role, does not contain any additional information, or that the conditional mutual information is strictly zero.

In the follow-up, priority should be given to improving the consistency of research questions, estimation objects and final conclusions, placing m6A as the main subject of research, multi-modification as exploration expansion, and perturbation as the source of subsequent hypotheses. The new analysis should be centered around a small number of pre-specified local patterns and verified using new independent data; if a reproducible pattern is not obtained, small effects and their applicable boundaries should be retained to avoid repeatedly cutting subgroups to look for positives.

## 9. Follow-up research directions

The next stage is recommended to focus on "conditional local structure correlation", and the main research questions are:> How stable is the association of RNA modifications with local structure in human publicly available transcriptome data; can the incremental information provided by structure be reproducible across conditions and technologies, controlling for sequence and measurement context?

Recent work can be divided into two steps. The first step is to organize the existing results into a complete research manuscript, taking the correlation, robustness, in vivo and in vitro differences, cross-technical review and prediction boundaries of m6A as the main lines, and the results of multiple modifications and perturbations as exploratory expansions. The second step prespecifies a small number of biologically meaningful stratifications, such as 5′UTRs, CDS, 3′UTRs, sequence motifs, and coverage levels, studies local structure curves and formal interaction effects, and double-checks on new independent data.

Perturbation studies should continue after the replicate layer data have been recovered and correct experimental unit inferences established. Ideally new data should provide both site-by-site modification changes and structural changes in the same sample and include rescue, inactive drug, or other controls that identify off-target effects. Once these conditions are met, the project is suitable to advance from Directionality to Causal evidence.

Before obtaining stable multi-modification and multi-condition corpus, it is not recommended to prioritize training Transformer, carry out large-scale parameter adjustment, or build large-scale interactive websites. Current limitations arise from data comparability and independent validation rather than insufficient model complexity.

## 10. Stage conclusion

ModStruct has completed four phases of research from single modification correlation to method robustness, cross-condition review, multi-modification exploration and perturbation directionality analysis. The existing evidence supports "small, condition-sensitive observational correlations between m6A and local experimental structures", but does not support stronger conclusions such as "experimental structures can stably improve modification predictions", "different modifications have formed clear structural niches" or "m6A changes have been shown to lead to structural remodeling".

The four phases of work can be concluded as a complete phase. The next wave of research should be less lateral, focus on finding condition-specific local structural patterns that can be replicated in independent data, and put experimental units, measurement techniques, and causal identification at the center of analytical design.

## Data, Code and Ethics Notes

This study used publicly available data and did not include newly recruited human participants or newly collected personal sensitive information from this project. Code, frozen configurations, data provenance records, machine-readable results, and stage reports are kept in this project directory. Use of external data remains subject to the original database and data provider license terms.

## Main evidence documents within the project

- [Original project idea](../../research_notes/chat_idea_organized.md)- [Fourth Phase Restart Scope List](../../research_notes/ModStruct_Restart Scope List.md)
- [GLORI field semantic verification](../../../metadata/provenance/01_glori_field_semantics.md)
- [Data source and file integrity verification](../../../metadata/provenance/source_verification.md)
- [HEK293T coordinates and initial intersection report](../../../results/reports/05_hek293t_initial_intersection_report.md)
- [HEK293T analysis data report](../../../results/reports/06_hek293t_analysis_dataset_report.md)
- [HEK293T Association Validation Report](../../../results/reports/07_hek293t_association_validation_report.md)
- [HEK293T prediction validation report](../../../results/reports/08_hek293t_prediction_validation_report.md)
- [HeLa Analysis Data Report](../../../results/reports/09_hela_analysis_dataset_report.md)
- [HeLa Association Review Report](../../../results/reports/09b_hela_association_report.md)
- [HeLa frozen model transfer report](../../../results/reports/09c_hela_model_transfer_report.md)
- [Phase 1 Method Expansion Summary](../../../results/reports/phase1_summary_report.md)- [Phase 2 Completion Report](../../../results/reports/22_phase2_completion_report.md)
- [Phase 3 Completion Report](../../../results/reports/27_phase3_completion_report.md)
- [Phase 4 Completion Report](../../../results/reports/31_phase4_completion_report.md)
- [Phase 4 Directionality Analysis](../../../results/reports/29_phase4_directionality_report.md)
- [Phase 4 Pre-registration and Cause and Effect Diagram Contract](../../../metadata/decisions/phase4_preregistration.md)

## External data and literature

- NCBI Gene Expression Omnibus. GSE264642: *DHX36 binding induces RNA structurome remodeling and regulates RNA abundance via m6A/YTHDF1*.
- Zhang, Y., Zhao, J., Chen, X., et al. (2024). DHX36 binding induces RNA structurome remodeling and regulates RNA abundance via m6A reader YTHDF1. *Nature Communications, 15*, 9890. https://doi.org/10.1038/s41467-024-54000-y