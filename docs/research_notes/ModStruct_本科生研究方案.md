# ModStruct: Correlation and predictive analysis of human m6A modification levels and local RNA structure

**Research positioning. ** This project is based on publicly available data to study the relationship between m6A modification levels and local RNA structure detection signals in human cells, and to evaluate the additional predictive value of experimental structural information relative to sequences and their computational structural features. Using bioinformatics, statistical modeling and lightweight machine learning methods, it is suitable for undergraduates with Python basics to complete under the guidance of instructors. Recommended cycle time is 12 weeks, with a 2-week buffer; week 8 forms the smallest version that can be delivered independently.

**The boundary between evidence and plan. ** Literature and data record verification as of September 9, 2026. What this article presents is the research design, not the results of the data analysis that has been run. Core papers, GEO projects, and some sample records have been checked; field-by-field inspection of data files, coordinate unification, or actual intersection calculations have not been completed, so there is no commitment to effective sample numbers, prediction improvements, or statistical significance. The screening thresholds, workload and acceptance thresholds below are all recommended settings for this project and must be frozen after the small-scale data audit and before the inspection results are officially viewed.

**1. Research questions and scope**

The recommended title is "Study on the association between m6A modification levels and RNA local structure based on open omics data". ModStruct As a project abbreviation, the English working title is available *ModStruct: Associations between m6A methylation levels and local RNA structural profiles*.

The main question is: **Does local RNA structure probing signal correlate with m6A modification ratios in human mRNA sites with reliable measurements; does this correlation persist after controlling for available sequence, transcript region, and measurement quality information and can it be supported in another cell line? **

Break the problem into three parts that can be delivered independently:

| Problem | Main Input | Main Output | Priority |
|---|---|---|---|
| Q1: Is the modification ratio of the detected m6A sites related to the neighborhood structure signal? | GLORI quantification, icSHAPE reactivity, background covariates | Effect size, confidence intervals, local curves | Must do |
| Q2: Does the experimental structure improve predictions of modification ratios? | Sequence, sequence predicted structure, experimental reactivity | Ablation comparison on the same test set | Must do |
| Q3: Are the above results robust across data sets? | HEK293T discovery set, HeLa review set | Correlation direction, model gain and difference explanation | Full version must do |
| Q4: Can the structure distinguish between modified and low-modified candidate sites? | Reliable positives and low-modification background with detection capabilities | Classification model, information gain | Only do it with qualified background data |

Q1 and Q2 default to "m6A sites that have been detected and passed quality control", and the response variable is a continuous modification ratio. This allows the project to not rely on treating unreported sites as negative. The cost is that the conclusions are only applicable to this selected group; detection thresholds, coverage, and selection bias still need to be discussed, and cannot be generalized to the rules of all A sites in the entire transcriptome.

The main lines for undergraduate students are limited to: human, mRNA, m6A, one discovery cell line and one review cell line, and published processed data. lncRNA, multiple modification comparisons, causal inference, basic model training, and database websites are not considered as closing conditions. The detection technologies and targets of different modifications are quite different, and initial merging will amplify the workload and make interpretation difficult.

**2. What has been demonstrated by existing research and what else can the project do**

This project is not starting from a blank slate. The icSHAPE work in 2015 has compared the structural characteristics of modified and unmodified sites in the same sequence motif and analyzed the changes after Mettl3 deletion in mouse cells; the m6A-switch study in the same year provided an example of the mechanism of local structural changes binding to HNRNPC. Therefore, "m6A is related to RNA structure" cannot be regarded as a new discovery in this project. [^1][^2]

In 2024, Shachar et al. combined hybridization systems, reporting experiments, and genetic variation analysis to study how sequence and structure affect m6A deposition. This also shows that both "structure affects modification" and "modification affects structure" have research foundations, and observational integration cannot determine the direction on its own. [^3]

Special attention also needs to be paid to SMART-m6A published on August 14, 2026. The study has combined sequence and computational structures to perform deep learning, ablation, quantitative prediction and cell line correlation analysis, and used icSHAPE to examine structural features. Its main model structure input comes from sequence prediction. Therefore, just adding a "sequence + structure" model or drawing a cell line difference curve are not enough to claim the originality of the method. [^4]

StructRMDB already provides resources for computing structural changes related to modifications. Therefore, mapping modifications to predicted structures and making them into web pages is not a priority for this project. [^5]

**The feasible contribution of this project is an empirical analysis with clear scope and reproducibility: with experimental reactivity as the core, sequence representation gain and experimental measurement gain are separately investigated, and the conditions for establishing the quantitative conclusion are quantified under strict division, coverage control and review with another data set. ** This is an appropriate undergraduate research objective; fulfillment of independent thesis requirements will depend on actual results and further comparison with existing work. Expressions such as "first time", "full disclosure" or guaranteed publication are not allowed.

Three equally qualified results can be formed: there is a stable gain in the experimental structure; the rough correlation is significantly weakened after controlling the background; and the gain only exists in specific measurement conditions or areas. Project quality depends on whether the evidence is reliable, not on whether the expected direction is obtained.

**3. Data selection and specific entrance**

| Purpose | Data | Recommended range | Verified information | Pending implementation |
|---|---|---|---|---|
| Discovery Set: Structure | GSE74353 | in vivo icSHAPE for HEK293T; in vitro is optional | GEO lists standalone InVivo, InVitro BaseReactivities compressed files | Transcript version, missing value encoding, whether to provide duplicate level scores |
| Discovery set: m6A | GSE210563 | HEK293T GLORI+ rep1/rep2 | Sample GSM6432590, GSM6432591; provide processed CSV | Column definition, complete measurable background, repeated integration method |
| Review set: Structure | GSE145805 | HeLa icSHAPE | Sample HeLaD/HeLaN; processed TXT archive | Processed reactivity must be used, D/N cannot be treated as biological duplicates |
| Review set: m6A | GSE210563 | HeLa hypoxia- GLORI+ rep1/rep2 | GSM6432595, GSM6432596; hypoxia- is the normoxic control for this experiment | Differences in culture conditions and structural research, whether there is a more suitable control under the same conditions |
| Cross-technology extension | GSE162356, part of GSE162357 | HeLa poly(A) m6A-SAC-seq | GEO lists processed BEDs for HeLa, HEK293, HepG2 | Availability of quantification fields, quality thresholds, and background sites |
| Search for alternative structural data | RASP v2.0 | Measured data consistent with the above cell conditions | Includes structure detection data, and includes some AI filling signals | Distinguish between measured and filled when downloading, trace back to the original research |

The data record entries are [GSE74353](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74353), [GSE210563](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210563), [ GSE145805](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145805), [GSE162356](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162356). These are candidate pairwise relationships, not common measurements of the same batch of cells. [^6][^7][^8][^9]

The first batch of recommendation documents are:

- `GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities.txt.gz`, the GEO marker is approximately 8.8 MB.
- `GSE74353_HS_293T_icSHAPE_InVitro_BaseReactivities.txt.gz`, approximately 7.0 MB, used only when expanding.
- GLORI processed files for GSM6432590, GSM6432591; first sample listed`GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv.gz`, about 7.7 MB.
- Download the required files for GSE145805 during the HeLa review stage. entire`GSE145805_RAW.tar`Marked at approximately 424.4 MB in GEO, the RAW in the name does not imply that it is all FASTQ; the record states that it is a TXT archive.

The file size is only used for download planning and is subject to the real-time directory. Not allowed`totalm6A.FDR`The file name presumes that it contains all tested negative sites. [^6][^8][^10]

**Pairing rules. **HEK293 and HEK293T are processed separately; normoxia and hypoxia, in vivo and in vitro are processed separately; different control groups within a study are not automatically combined. Hypoxia- of HeLa is a control condition for a particular study and cannot be written as "exactly the same conditions" as another experiment. The same cell line can only reduce part of the differences, but cannot eliminate differences in batches, culture, passages, and isoform composition.

Human RNA MaP can be used as a follow-up resource, but does not enter the first stage: it mainly corresponds to U2OS and does not match the modified samples prioritized in this protocol. Their DMS-TRAM-seq work remains marked as preprint in the PubMed record for this review and should be cited as such. [^11]

An early data note must also be corrected: GSE50676 corresponds to Wan et al.'s human family RNA structure study using PARS and cannot be classified as DMS-seq; this plan does not mix it into the icSHAPE main analysis. [^12]

**4. Data audit and minimum viability check**

The goal of weeks 1-2 is to verify whether a reliable site-level union table can be built, rather than to pursue model scores. First read the header and a small number of records of each file, and then parse the entire table. Record the source link, file name, download date, SHA-256, cell line, processing conditions, experimental method, reference version, coordinate rules and column descriptions.

Each type of document must answer the following questions:

1. Does a row represent a genomic site, a transcript site, an interval, or a complete transcript?
2. Do the coordinates start from 0? Is the interval closed on the left and open on the right? How is chain information represented?
3. Are the structure values ​​normalized reactivity, raw counts, predicted values, or AI imputed values?
4. Is the m6A value modification proportion, mutation proportion, significance or confidence score? Is it background corrected for chemical conversion?
5. Which rows have been filtered for significance? Is there a reliable measurable background for all A sites?
6. Are biological duplicates still separable, or have the authors been merged?

The suggested feasibility judgment is as follows. They are engineering thresholds, not statistical efficacy guarantees:

| Audit results | Decision-making |
|---|---|
| The discovery set has about 2,000 available sites, involving about 300 genes, and the distribution and quality are basically reasonable | Advancing the complete main line |
| About 500-2,000 sites or less gene coverage | Reduce the number of features, focusing on regression and description; reduce the tree model to optional |
| Less than about 500 sites, or some strata have few independent genes | Check mapping and coverage first; if necessary, change to a pre-listed data combination |
| Unable to determine coordinate or quantitative field meaning | Suspend formal statistical analysis of this data without guesswork |
| No detection coverage information for low-modification background | Keep the main line of continuous ratio and cancel the modified/unmodified classification |

A screening flow chart is generated before the end of the second week: number of downloaded records → site deduplication → reference sequence verification → mRNA annotation → structural coverage → repeat QC → final number of sites and genes. The actual calculated value must be reported for the sample size, and the overall database size cannot be used instead of the intersection size.

**5. Coordinate unification and joint data table**

Priority is given to GRCh38 used by GLORI records, but the conversion relationship must be established based on the transcript annotation actually used in the structure file. Don't just make a direct connection based on "all are human" or "all are hg38". The transcript version of GENCODE/Ensembl affects the splicing path and position, so you cannot just remove the version number and merge blindly. [^10]

The processing sequence is: check the original transcript sequence and annotation of the structure file → map its position to the stranded genomic coordinates → check against the GLORI genomic site → extract the spliced ​​local RNA sequence by the selected transcript. If cross-version migration is a necessary step, re-verify the central base, neighbor sequence and chain direction after migration. The genomic reference base for the reverse strand may be T, but the central base in the direction of the transcript must be A.

Taking the unique "genome site + chain + cellular conditions" as the core biological unit. When the same site is mapped to multiple isomers, they should not be treated as multiple independent samples. Records with unambiguous local sequence and region assignments were retained in the main analysis; the remainder were counted separately and included in sensitivity analyzes or excluded. If the representative transcript rule is adopted, the rule must first be fixed and explained that the problem of short-read isoform mixing cannot be solved. Isomer perception analysis is not required for undergraduate students.

Suggested union table fields:

| Field Group | Content |
|---|---|
| identification | site_id, gene_id, transcript_id, chr, position, strand, cell_line |
| source | structure_accession, modification_accession, replicate, condition, reference_version |
| response | m6A_fraction, raw molecule/read count (if any), coverage, raw confidence metric |
| Sequence | Center 201 nt sequence, DRACH subtype, GC ratio, 3-mer frequency |
| Comments | CDS/5′UTR/3′UTR, position from stop codon and splicing boundary, expression level (if any) |
| Experimental structure | -50 to +50 nt reactivity array, rms mask, individual window coverage |
| Predicted structure | Pairing probability, MFE, pairing state entropy, etc. of 201 nt sequence prediction |
| Analysis flags | Reasons for inclusion/exclusion, missing field flags, gene_group, split_id |

At least 20 sites were manually checked, including positive strand, negative strand, exon boundaries and multiple isomers. Formal scripts need to have coordinate assertions and length checks that run automatically. The inspection here is a necessary verification to prevent bioinformatics misalignment, rather than proving the overall biological conclusion with a few manually selected sites.

**6. Quality control, response variables and missing values**

GLORI provides single-base quantitative m6A information, but ratio estimates are still affected by effective coverage and chemical conversion background. Prioritize the quality control defined in the original text and file; only add coverage filtering when the corresponding fields are provided in the data. The ordinary mutation rate must not be renamed uncorrected to the modification rate. [^13]

The main analysis prioritized retaining loci whose two replicates met the author's quality requirements. Inter-replicate correlations, absolute difference distributions, and coverage relationships are examined separately; not just the correlation coefficient is reported. If there is only a proportion, the replicate mean is used by default and each replicate result is retained for sensitivity analysis; if there are comparable effective counts, the counts can be combined according to the original process, but the technical reads cannot be treated as independent biological replicates.

If you have a clear valid depth field, you can compare the retention rates of ≥10, ≥20, and ≥50 during audit, and select the main threshold before looking at the formal correlation results. Don’t pick a coverage threshold to get significant results. When the modification table does not have a depth field, it means that it can only inherit the author filtering and cannot claim to have controlled the depth of site detection.

The main window of the experimental structure is defined as 10 nt on each side of the center, excluding site 0; the main index is calculated only when the effective reactivity on both sides is at least 70%. The number of valid samples for each relative position was recorded individually for the ±50 nt window used for presentation. all`NULL`、`-999`The missing codes are converted to missing according to the original document instructions and are not filled in with zeros.

The main analysis retains the measured structure values; AI padding signals are not mixed with the measured results. RASP v2.0 does provide deep learning imputation for some data, so this distinction must be preserved after downloading. [^14] When covariates are missing, lightweight imputation is only fitted to the training set; coverage related to missing structures should be used as a quality diagnostic. Just because a model can handle NaN does not mean that missingness bias has disappeared.

**7. Structural features: few but clear**

Structure detection reactivity is a structure-related measure, not a direct pairing probability. Reports are collectively referred to as "icSHAPE reactivity" or "experimental structural signal"; measurement boundaries are noted when interpreted as evidence of local flexibility or accessibility. Probe accessibility is also related to the local chemical environment and molecular contacts, and high reactivity cannot always be explained as a specific stem being opened. [^1][^15]

| Characteristics | Definition | Role |
|---|---|---|
| R_flank10 | Mean of -10…-1 and +1…+10 effective reactivity | Q1 only major structural indicator |
| R_up20, R_down20 | 20 nt mean upstream and downstream, excluding the center | Explore neighborhood directionality |
| R_far | The mean of -50…-21 and +21…+50 | Distinguish between nearest neighbors and further background |
| R_0 | Central site reactivity | Additional analysis only to avoid mainline reliance on central readout |
| Coverage | Proportion of valid values ​​in each window | QC and deletion bias diagnosis, not packaged into biological structural features |
| P_unpaired | Calculate the predicted central or neighborhood unpaired probability | Sequence prediction structure comparison |
| MFE_per_nt | Fixed window MFE divided by length | Compute features, not interpreted as modified energy effects |
| H_pair | Compute the pairing state uncertainty of the ensemble | Information theory extension |

The predicted structure uses a fixed 201 nt sequence window, consistent with the main model input sequence range; the first version uses only standard A/C/G/U parameters. Use ViennaRNA's partition function to calculate pairing probabilities instead of deriving probabilities from just a dot-bracket structure. RNAfold`-p`This type of calculation is supported. [^16]

If p(i,j) is defined as the predicted probability of pairing site i with j, u(i)=1−Σj p(i,j), then it can be defined:

$$H_i=-u_i\log_2u_i-\sum_jp(i,j)\log_2p(i,j).$$

It is agreed that 0·log(0)=0. These states are mutually exclusive; overlapping labels such as paired, unpaired, stem, and loop cannot be placed in the same probability distribution to calculate entropy. H_i is the pairing state uncertainty for a given model, window, and temperature and is not a direct measure of the conformational population or dynamics within the cell.

If the preparation of the folding software environment takes more than half a day, complete the experimental signals and sequence baselines first, and then add the calculation structure. Do not modify core dependencies in the existing PyTorch environment for this project. Window truncation may ignore long-range pairings and should be stated in the Discussion; window sensitivity analysis of 401 nt was only done when available and simultaneously extended sequence baselines.

**8. Statistical analysis: answer the correlation first, then explain the model**

**8.1 Descriptive diagram. ** Plot modification proportion distribution, structure coverage distribution, repeat consistency, and ±50 nt average reactivity curves grouped by modification proportion in the discovery set. Grouping boundaries are determined and frozen in the discovery set; continuous scale remains the official response, and grouping is used for visualization only. The number of effective sites is provided below the curve to avoid mistaking coverage changes as structural changes.

**8.2 Primary Relevance. ** Let y be the proportion of m6A after repeated integration and R be R_flank10. Fit a low-dimensional model:

$$y_i=\beta_0+\beta_RR_i+\gamma^TC_i+\epsilon_i.$$

C pre-includes the DRACH subtype, transcript region, local GC, position from the stop codon, as well as obtainable expression levels and coverage metrics for both experiments. Primary analyzes restricted to conventional DRACH must specify their scope; non-DRACH sites may be described separately. It cannot be deduced from deletion of rare subtypes or atypical sites that they have no biological significance.

Mainly report β_R, modification proportion change corresponding to one standard deviation of structural signal change, 95% confidence interval and sample/gene number. If using a linear model, check the residuals and apparent nonlinearities; it is an average correlation approximation, not a modified probability generation mechanism. For significantly curved relationships, a low-degree-of-freedom spline can be preset for sensitivity analysis instead of continuing to increase the degrees of freedom until significant.

Main confidence intervals were bootstrapped clustered by gene, suggested 1,000 times, and compared to clustered robust standard error results. All sites of a gene are entered together in the same resampling; the effective sample size cannot be inflated with independent site bootstrap. The results with equal weight per gene were added as a sensitivity analysis to check whether a few long transcripts dominate the conclusions.

**8.3 Multiple comparisons. ** The main test of R_flank10 is fixed in advance; other structural windows, regional hierarchical and position-by-position tests define test families respectively, using Benjamini–Hochberg FDR. If the position-by-position plot only has pointwise intervals, make sure it is not the simultaneous confidence band of the entire curve. A string of adjacent salient locations cannot be treated as many independent discoveries.

**8.4 Sensitivity analysis must be done. ** Compare at least: exclusion of center 0 vs. further exclusion of -2…+2; only high structure coverage records; use of two GLORI repeats each; equal weighting by gene; restricted to sites without isoform ambiguity. When the sample is sufficient, add CDS and 3′UTR layering. If the correlation disappears after removing the signal near the center, the conclusion is narrowed to local detection features rather than proving the existence of large-scale refolding. [^15]

**8.5 Unmeasured factors. **Protein binding, transcription status, cell culture differences, and isoform mixing may not necessarily be adequately controlled. "Adjusted association" in a regression simply means adjusting for the listed covariates. Don't write "all confounding is ruled out" or "prove that structure independently determines modification."

**9. AI module: Use ablation to determine where information comes from**

The main task of the AI ​​is to predict the continuous m6A modification ratio. Only two types of models are compared: Ridge regression as a stable baseline and HistGradientBoostingRegressor as a nonlinear model. The latter is an implementation of gradient boosted trees, suitable for tabular features; the mainline does not require PyTorch, Transformer, or GPUs. [^17]

Each model uses the same response, the same batch of sites, and the same data partitioning. Define X as the sequence feature of 201 nt (one-hot plus preset 3-mer frequency), P as the predicted structure of the same sequence, R as the above experimental reactivity feature, and C as the common background covariate.

| Model Number | Input | Answered Questions |
|---|---|---|
| M0 | C | How much can the regional and measurement context itself explain? |
| M1 | C+X | How strong is the sequence baseline? |
| M2 | C＋X＋P | Does computing structural representation help this learner? |
| M3 | C＋X＋R | What is the additional predictive value of experimental structures relative to sequences? |
| M4 | C＋X＋P＋R | After the calculation structure is already available, can the experimental structure still provide gains? |

**The main model comparison is M4 vs. M2, and the main metric is the MAE difference on the independent test set: **

$$\Delta MAE=MAE(M2)-MAE(M4).$$

Positive values ​​indicate that the error decreases after adding the experimental structure. Secondary measures are Spearman correlation, RMSE, and R²; negative R² is allowed. If the ratio is expressed as 0-1, the MAE multiplied by 100 corresponds to the percentage, and the unit must be written when reporting. Don't use "accuracy" to describe continuous proportional regression.

Ridge recommends comparing only a small number of regularization strengths, such as 0.1, 1, 10, and 100; the lifting tree presets about 6-8 sets of parameters and does not conduct large-scale searches. The default random internal early stopping partitioning of tree models may span genes and needs to be disabled or an explicit group validation scheme used. All scaling, interpolation, and parameter selection occur only in development data.

When interpreting the model, first check whether it is more effective than the baseline before doing feature group ablation or test set permutation importance. For highly related structural windows, give priority to the whole set of explanations to avoid fabricating independent mechanisms for each similar feature. Replacement will produce unnatural sequence-structure combinations, so the main evidence is still model comparison after retraining; the importance of replacement is only used for model diagnosis. SHAP is not a must-do, nor is it a cause-and-effect explanation. [^18]

**10. Data division and external review**

The main partition is in genes: ~80% of HEK293T's genes as a development set and 20% as a one-time internal test set. Use 50% off GroupKFold tuning in the development set. If the number of effective genes is small, reduce it to 3 fold and reduce the model size. The purpose of GroupKFold is to prevent a group from entering training and validation at the same time. [^19]

The division also deals with cross-gene repeat sequences: identical 201 nt windows must be in the same group, or only one copy is retained; high similarity windows are audited and grouped together by sequence clusters and gene relationships when possible. Grouping by genes alone does not mean completely eliminating homologous sequence leakage. All repetitions of match sets and sites must also remain in the same group.

The master model reports final results once on an internal test set; test performance cannot in turn be used to select windows, features, or thresholds. The ΔMAE interval was calculated using paired bootstrap with test genes as units, and all models sampled the same batch of genes each time. This interval mainly reflects the sampling uncertainty of the fixed model on the test genes, and does not include all experimental batches and training randomness. You can use a small number of preset seeds to supplement training stability, but you can't pick the best ones to run.

HeLa review is divided into two levels:

- **Association Review:** Effects are estimated independently in HeLa using frozen R_flank10, covariates and QC rules. Report the direction, magnitude, and interval; if it is contrary to the discovery set or the interval is wide, write it directly.
- **Migration evaluation:** Freeze the HEK293T model and preprocessing, and test directly in HeLa without HeLa label parameter adjustment. Report all qualified sites separately, as well as a strict subset after excluding genes/repetitive sequences that appeared in the development stage; if the latter is insufficient, it can only be called cross-cell background migration, not generalization of new genes.

The experimental years and processing procedures of GSE74353 and GSE145805 are different. Therefore, cross-cell line variation and study batch variation are entangled here; HeLa review failure cannot be attributed to "cell specificity" alone, nor does success prove to be universal for all cells. When necessary, only effect directions and relative model gains within each data set were compared, without merging raw reactivity. Batch calibration determined by test labels is prohibited.

**11. Information Theory Module: Optional, but must be defined accurately**

The "information" in this project can have three specific meanings: the integration of sequence and experimental data; the calculation of the pairing state entropy of the ensemble; and the degree of reduction in prediction uncertainty after adding structure to the model. The first two items can be entered directly into the main plot, and the third item serves as an expansion for weeks 9-10.

There is a key mathematical bound: if the predicted structure P is completely determined by the complete input sequence X and a fixed algorithm, that is, P=f(X), then theoretically:

$$I(Y;P\mid X)=0.$$

This does not prevent M2 from being preferred over M1, since limited data and limited models may be easier to use structural representations. Such an improvement should be called a "structural representation gain" and cannot be interpreted as P bringing new observational information beyond the complete sequence. If the model only uses compressed k-mer features, the predicted structure may also compensate for the sequence information lost during compression, so the results cannot be regarded as proof of independent biological information. The experimental reactivity R is not calculated solely from the input sequence and deserves separate evaluation.

It is not recommended for undergraduate students to directly estimate high-dimensional I(Y;R|X) in the first version: high-dimensional conditional independence testing requires clear assumptions, and the estimation results are easily affected by model errors and sparsity. [^20]

An easier-to-implement option is to define low, medium, and high modification levels based on the frozen proportion quantiles of the discovery set within the detected sites, and use classification as an auxiliary task; or perform secondary classification when qualified low modification background exists. Convert the reduction in model logarithmic loss to bits:

$$\Delta_{bits}=\{LogLoss(M2)-LogLoss(M4)\}/\ln 2.$$

The quantity is calculated on independent test samples and intervals are given by gene. It is called "the logarithmic loss gain under this prediction task", do not directly name it as the accurately measured biological condition mutual information; only when the model fully approximates the true condition distribution and meets the corresponding sampling conditions, there is an ideal relationship with the conditional entropy difference. Grouping thresholds, class proportions, and selection bias all change it. If negative, report it as it is and cannot be truncated to zero.

**12. Optional modification/low modification control analysis**

This part is not a prerequisite for the main plot. Classification labels are only constructed if there is a "tested and sufficiently covered" A site background. A outside the list of public significant sites is uniformly regarded as unlabeled, and y=0 cannot be automatically set; the non-detection of a site does not mean that it is truly unmodified.

Qualified controls need to match: central base A, DRACH subtype, transcript region, sequence GC, expression and two types of experimental coverage. Prioritize isogenic or isogenic backgrounds; positives without suitable controls are excluded and reported. The low modification threshold is determined based on the technical detection limit and proportion interval. A line with no counts cannot be arbitrarily called "less than 5%".

If 1-3 matching backgrounds are used for each positive, it is to obtain research samples with comparable conditions, not natural population sampling. AUPRC, calibration, and log loss must be reported along with the test set positivity ratio; AUPRC on the 1:1 balanced set cannot be directly compared to the naturally sparse candidate set. The model learns the detection difference under the label and sampling rule, not the true modification probability of all RNA molecules.

When low-modification backgrounds are unavailable, the best alternative is to complete a continuous proportion analysis and its restriction specification, rather than ad hoc training of a high-scoring "positive versus random A" classifier.

**13. Result chart and delivery content**

It is recommended that 5 main pictures be formed in the end; each picture answers a clear question, and no hypothetical experimental results are drawn in advance.

| Figure | Content | Questions to answer |
|---|---|---|
| Figure 1 | Data source, sample conditions, screening process, final number of sites/genes | How does the data come from, and how much is available? |
| Figure 2 | Repeat consistency, coverage distribution, local curves for different modification levels | Is the signal reliable and what is the rough correlation? |
| Figure 3 | Forest plot of adjusted effects and main sensitivity analysis | Does the association persist after controlling for background? |
| Figure 4 | The MAE and ΔMAE intervals of M0—M4 in the fixed test set | What is the actual contribution of the experimental structure? |
| Figure 5 | HeLa review, strict and non-strict migration subsets, representative sites | Can the conclusion be reviewed, and where does it fail? |

The representative site selection rules are first fixed, such as displaying 3–5 test sites from high coverage, no mapping ambiguity, including at least one model failure case. Don’t just pick the best-looking genes. If the predicted secondary structure is shown in the figure, it is clearly marked as a computational prediction and is not directly analyzed as an experiment.

The final delivery package includes: research report, the above main picture, source and data dictionary, filter count table, fixed partition table, analysis script, configuration file, environment version, operation instructions and restriction list. Negative or unstable results are also retained in the report. No need to make an online database.

**Fourteenth and 12th Week Arrangement and Acceptance Criteria**

A student with basic knowledge will spend about 12-15 hours per week, and the tutor will discuss planning for about half an hour per week. For beginners or those with busy courses, a 14-16 week schedule is more secure; the following is an estimate of workload, not a guarantee of completion time.

| Weeks | Work | Result that must be delivered |
|---|---|---|
| Week 1 | Read core papers, organize data and fields, and download processed samples | Literature comparison table, data list, file header examples |
| Week 2 | Coordinate mapping, small-scale intersections, coverage audit | Site inspection records, preliminary screening form, OK to continue/narrow down |
| Week 3 | Full union table, duplicate QC, deletion handling | Freeze data dictionary, determine number of available sites and genes |
| Week 4 | Local features, Figures 1-2, freeze major issues and divisions | Written analysis and configuration; thereafter the main plan will not be changed based on test results |
| Week 5 | Association Models, Gene Cluster Intervals, Sensitivity Analysis | Figure 3 and Effect Scale |
| Week 6 | Sequence Baseline, Predicted Structural Characteristics, M0—M4 Ridge | Reproducible Development Set Comparison |
| Week 7 | Lightweight tree model, melting, freezing model | Final parameters and model description |
| Week 8 | One-time internal testing, finishing the minimum version | Figure 4, complete single cell line report, runnable code |
| Week 9 | HeLa review data integration and same rule analysis | Review QC, correlation effect table |
| Week 10 | External migration, optional information gain, or cross-technology validation | Figure 5; Select only one extension |
| Week 11 | Organize explanations, failures and limitations | First draft of report, all illustrations and sources |
| Week 12 | Rerun key analysis from fixed intermediate table, instructor verification | Final report, code, defense material content outline |

**Week 8 Minimum Closing Criteria:** Confident joint data from one cell line; clear association analysis; fixed test comparison of sequence to experimental structure; coordinate and partitioning checks; sufficiently clear statement of limitations. When the predicted structure is not completed on time, M3 versus M1 may be compared as a minimal version, but no claim shall be made that computational structural and experimental structural contributions have been separated.

**Full Version Standard:** Completion of M4 versus M2 and HeLa review on Minimum Version. Do not use p<0.05, ΔMAE>0, or a certain AUC threshold as closing requirements. Insufficient samples, missing gains, or external failures can form a complete undergraduate research report as long as there is sufficient analysis and explanation.

**15. Computing resources and engineering organization**

The main process uses the processed form as input. It is recommended to use a computer with 16 GB memory to start and 32 GB for a more comfortable computer, and reserve about 10-20 GB project space; these are conservative plans and should be based on actual files and intermediate results. The mainline does not require a GPU, does not rerun the entire set of FASTQ in batches, and does not build a self-built RNA basic model.

Python works with pandas/numpy, scipy/statsmodels, scikit-learn, and regular plotting libraries; structure folding prepares ViennaRNA separately. The project's known callable Python is`E:\ancd\envs\my_pytorch\python.exe`, but whether the specific library is available should still be subject to environmental inspection. These dependencies have not been installed or modified in this program preparation.

The recommended directory is as follows; this is a future implementation structure and does not mean that these data and codes have been generated:

```text
ModStruct/
  README.md
  config/analysis.yaml
  metadata/manifests/source/source_manifest.csv
  metadata/provenance/data_dictionary.md
  data/raw_processed/
  data/interim/
  data/final/site_table.parquet
  results/tables/08_hek293t_group_assignments.csv
  src/03_audit_inputs.py
  src/04_validate_coordinate_mapping.py
  src/03_build_features.py
  src/04_association.py
  src/05_models.py
  src/06_validate.py
  results/tables/
  results/figures/
  docs/reports/final/
```

The configuration is at least fixed: data version, inclusion criteria, missing reactivity encoding, window, primary response, QC threshold, random seed, gene grouping, model parameter candidates, primary endpoint, exploratory test family. Notebooks are used for exploration, and formal output is generated by scripts reading configurations. Students should be able to explain the origin of each column and how each graph was calculated.

**16. Risks, Alternative Routes and Conclusion Boundaries**

| Risk | How to find out | How to deal with it |
|---|---|---|
| Coordinate or version mismatch | The center is not A, the motif is abnormal, and the intersection is extremely low | Return to the original sequence and annotation; do not solve it by adjusting the statistical threshold |
| The positive table has no measurable background | The file only lists significant sites and does not have a complete coverage table | Continuous proportion main line; cancel negative classification |
| Structural coverage is heavily biased toward high expression | The inclusion rate changes significantly with expression/coverage | High coverage sensitivity, indicate the applicable group, and narrow the scope if necessary |
| The signal is concentrated only in the modified ortho positions | The effect disappears after excluding -2…+2 | Described as a local detection feature, no global structural remodeling is inferred |
| Cross-training test of the same gene | split audit failed | Re-divided into groups, discarding leaked version performance |
| HeLa review failed | The direction is inconsistent or the error gain disappears | The reported batch cannot be separated from the biological conditions and is not selectively deleted |
| The calculated structure is effective and the experimental structure has no gain | M2>M1, but M4 is not better than M2 | It is expressed as the representation gain of the current baseline, which cannot prove that the experimental structure is useless |
| The range of the proportion within the positive is too narrow | The distribution is obviously cut off, and the repetitive noise is close to the site difference | Narrow the conclusion; do not arbitrarily filter the high/low end to improve predictability |
| The workload exceeds the semester | The joint table cannot be formed in the 4th week | Retain a cell line and measured reactivity, delete the tree model/information theory extension |

All results should be limited to the cells, measurement techniques, site filters, and sequences studied. Multiple genes provide statistical information, but do not equate to multiple independent cell experiments; genetic bootstrap is not a substitute for biological replicates. Neither classification nor regression alone can determine whether "structure affects modification" or "modification affects structure."

For future mechanistic validation, local structures could be compared with modified/unmodified RNA of the same sequence in a separate project, or by combining modified enzyme perturbation, rescue, and structural measurements. Enzyme perturbation itself may also change expression and cell status, and corresponding controls are still required. The undergraduate version does not rely on these labs for completion. [^2][^3]

**17. Project difficulties and work boundaries assisted by AI**

Most of the difficulties in this project are not in programming, but in data auditing, judgment and analysis discipline. They are listed below in descending order of error risk: the first two can only be solved by actually reading files and manual judgment, the middle two rely on rule freezing and conscious compliance, and the last defines how AI tools are used in the project.

**(1) Data alignment and coordinate unification (highest risk). ** GLORI and icSHAPE come from different studies, different years, and different transcript annotation versions. When they are combined into the same "site", any error in the coordinate starting method, interval opening and closing, genome version, transcript version, chain direction, and isoform attribution will invalidate all subsequent analysis as a whole, and the error is silent - the code is run, the picture is shown, but the conclusion is wrong. The manual verification and coordinate assertion required in Section 5 are not optional optimizations, but necessary hurdles to prevent "everything is wrong".

**(2) The true status of the data is unknown. ** The file has not been checked field by field when the protocol was prepared. Questions such as whether a row represents a site or the entire transcript, whether a column represents a modification ratio or a mutation rate, and whether duplicates have been merged by the author have not yet been answered. These can only be solved through real downloads, header readings, and column-by-column audits, and cannot be replaced by speculation or "AI-generated understanding." Any fields that have not been personally verified shall not be deemed to have been verified.

**(3) Prevent "adjusting for good-looking results". ** Choosing coverage thresholds, choosing the best from multiple random seeds, allowing the same gene to be tested across trainings, and only displaying good-looking sites are all analytical discipline issues rather than technical issues. There are no code error prompts and can only be avoided by relying on the rules of freezing in Sections 4, 9, and 10 and consciously avoiding them.

**(4) Explain the boundaries. ** Reactivity is not equal to pairing probability, "adjusting covariates" is not equal to "excluding all confounding", and the gain when the predicted structure is completely calculated from the sequence can only be called "representation gain". These boundaries do not affect whether the diagram can be produced, but they determine whether the conclusion is tenable.

**(5) Core difficulty will not be reduced by AI tools. ** Generative AI can significantly lower the threshold for writing scripts, checking information, drawing, and writing documents, but it cannot complete three things for people: actually reading data files, making judgments about the meaning of fields, and maintaining analytical disciplines without being broken by the inertia of "making the results beautiful." What needs to be more vigilant is that AI tools will adapt to cues and continue to optimize until the results look good, allowing errors such as p-hacking and data leakage to be covered up by smooth report packaging. Therefore, this project positions AI as an assistant for writing code and checking information. The students and tutors are jointly responsible for the meaning of data, analysis and judgment, and the boundaries of conclusions. The use of AI will be disclosed in accordance with the school's requirements during the official opening and closing of the project.

**18. Proposal summary and research value**

The abstract used to develop the title can be rewritten directly:

> m6A modification of RNA correlates with local structure, but correlations in published data may be affected by both sequence composition, transcript region, and measurement coverage. This project intends to integrate the GLORI quantitative modification data of human HEK293T cells and the icSHAPE experimental structural data, and use the detected and reliable m6A sites as the object to study the relationship between modification ratio and local structural reactivity. The project will build traceable site-level joint data, employ background adjustment and statistical analysis by gene clustering, and assess additional predictive value through model ablation of sequences, predicted structures, and experimental structures. Further review the correlation and migration performance with HeLa data to clarify the dependence of the conclusion on the cell background, experimental batch and detection range. It is expected to deliver reproducible analysis procedures, structural correlation results, and undergraduate research reports to provide empirical evidence for distinguishing the value of RNA sequence characterization and experimental structural information.

The subject affiliation can be written as “interdisciplinary research on bioinformatics, statistical learning and RNA epitranscriptomics”. The role of AI is to create a more interpretable prediction model; the role of information theory is to define the predictive information provided by different inputs and its boundaries. This project does not require a commitment to develop new neural networks, nor does it require packaging general correlational analysis into causal mechanisms.

**19. Instructions for searching and document usage**

This protocol uses a targeted range search, covering m6A—structural mechanisms, quantitative modification data, RNA detection resources, structure-assisted prediction, and model validation methods; it is not a PRISMA systematic review, nor does it claim to be an exhaustive search. Main search terms include`m6A icSHAPE structure association`、`GLORI GSE210563`、`GSE74353 HEK293T`、`m6A structure prediction 2026`、`StructRMDB`、`RASP v2.0`and`conditional independence testing`。

Priority will be given to original research from journals, official GEO records, and official software documents. Preprints are individually annotated, and PMC inclusion is not considered proof of peer review. Some database pages failed to access intermittently. The core metadata has been verified through readable official records, but it cannot be claimed that the actual data content or download links have all been verified. Source sizes and technical descriptions do not replace actual measured intersections for this project.

This article was organized using AI-assisted search and protocol drafting; no experimental results, analytical performance, or statistical significance were generated or assumed to be actual results. When formally opening and closing the topic, students and tutors should read the core methods and data descriptions, and disclose the use of AI as required by the school.

**Sources and References**

[^1]: Spitale RC, Flynn RA, Zhang QC, et al. *Structural imprints in vivo decode RNA regulatory mechanisms*. Nature. 2015;519:486–490. [Original](https://pmc.ncbi.nlm.nih.gov/articles/PMC4376618/), [DOI](https://doi.org/10.1038/nature14263). This paper exists [2015 Corrigendum Record](https://pubmed.ncbi.nlm.nih.gov/26416736/); its charts should be checked together before formally reproducing them, and this plan will not reuse its specific values ​​as results.

[^2]: Liu N, et al. *N6-methyladenosine-dependent RNA structural switches regulate RNA–protein interactions*. Nature. 2015. [Original](https://www.nature.com/articles/nature14234). Used to illustrate existing evidence for a structural switch mechanism without generalizing this mechanism to all m6A sites.

[^3]: Shachar R, Dierks D, Garcia-Campos MA, et al. *Dissecting the sequence and structural determinants guiding m6A deposition and evolution via inter- and intra-species hybrids*. Genome Biology. 2024;25:48. [Original](https://pmc.ncbi.nlm.nih.gov/articles/PMC10870504/), [DOI](https://doi.org/10.1186/s13059-024-03182-1).

[^4]: Sun M, Zhang D, Li Z, Lin Y. *Structure-aware deep learning enhances m6A prediction and reveals cell type-associated RNA structural signatures*. PLOS Computational Biology. 2026;22(8):e1014649. 2026-08-14. [Original text](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014649), [Code](https://github.com/mzsun-ai/SMART-m6A).

[^5]: Zhang Z, Wang X, Zhou J, et al. *StructRMDB: A database of RNA modification sites that affect RNA secondary structure*. Computational and Structural Biotechnology Journal. 2025;27:5493–5502. [Publication Record](https://pubmed.ncbi.nlm.nih.gov/41439023/), [Original](https://pmc.ncbi.nlm.nih.gov/articles/PMC12720313/), DOI: 10.1016/j.csbj.2025.11.058.

[^6]: NCBI GEO. *RNA duplex map in living cells reveals higher order transcriptome structure*, GSE74353. [Data record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74353). This solution only selects the icSHAPE files of HEK293T and does not classify the entire project as icSHAPE.

[^7]: NCBI GEO. *Absolute quantification of single-base m6A methylation in the mammalian transcriptome*, GSE210563. [Data record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210563); [HeLa normoxic control sample GSM6432595](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6432595).

[^8]: NCBI GEO. *RNA secondary structome of different cell lines*, GSE145805. [Data record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145805). Relevant original research: *Predicting dynamic cellular protein–RNA interactions by deep learning using in vivo RNA structures*. [Full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC7900654/).

[^9]: NCBI GEO. *Quantitative and whole-transcriptome N6-Methyladenosine profiling at single-nucleotide resolution [polyA_RNA]*, GSE162356. [Data record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162356). Related research: Hu L, et al. *m6A RNA modifications are measured at single-base resolution across the mammalian transcriptome*. Nature Biotechnology. 2022;40:1210–1219. [Full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC9378555/), DOI: 10.1038/s41587-022-01243-z.

[^10]: NCBI GEO. *HEK293T_GLORI+_rep1*, GSM6432590. [Sample and processing documentation](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6432590).

[^11]: Farenhem K, Whitfield TW, Chouloute A, Jain A. *Transcriptome-wide RNA accessibility mapping reveals structured RNA elements and pervasive conformational rearrangements under stress*. bioRxiv. 2025. DOI:10.1101/2025.06.05.658101. [PubMed preprint record](https://pubmed.ncbi.nlm.nih.gov/40502122/), [Full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC12157405/), [Human RNA MaP](https://humanrnamap.wi.mit.edu/).

[^12]: Wan Y, Qu K, Zhang Q, et al. *Landscape and variation of RNA secondary structure across the human transcriptome*. Nature. 2014;505:706–709. [Full text and PARS Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC3973747/), [Publication Record](https://www.nature.com/articles/nature12946).

[^13]: Liu C, Sun H, Yi Y, et al. *Absolute quantification of single-base m6A methylation in the mammalian transcriptome using GLORI*. Nature Biotechnology. 2023;41:355–366; Published online on 2022-10-27. [Publication and Abstract Record](https://pubmed.ncbi.nlm.nih.gov/36302990/), DOI: 10.1038/s41587-022-01487-9.

[^14]: Mu K, Fei Y, Xu Y, Zhang QC. *RASP v2.0: an updated atlas for RNA structure probing data*. Nucleic Acids Research. 2025;53(D1):D211–D219; published online in 2024. [Original](https://academic.oup.com/nar/article/53/D1/D211/7901281), [Database](https://rasp2.zhanglab.net/).

[^15]: *RNA infrastructure profiling illuminates transcriptome structure in crowded spaces*. Cell Chemical Biology. 2024;31(12):2156–2167.e5. DOI:10.1016/j.chembiol.2024.09.009. [Original](https://www.sciencedirect.com/science/article/pii/S2451945624004057), [PubMed](https://pubmed.ncbi.nlm.nih.gov/39447577/). Boundaries used to define probe accessibility, local chemical environment, and interpretation of m6A neighbor signals.

[^16]: ViennaRNA official documentation. *The Program RNAfold*; *Predicting various Thermodynamic Properties*. [RNAfold Tutorial](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/tutorial/RNAfold.html), [Partition Functions and Probabilities](https://www.tbi.univie.ac.at/RNA/ViennaRNA/refman/partfunc/thermodynamics.html). Accessed 2026-09-09.

[^17]: scikit-learn official documentation. *Histogram-Based Gradient Boosting*. [Model Description](https://scikit-learn.org/stable/modules/ensemble.html). Access date 2026-09-09; local actual version recorded during implementation.

[^18]: scikit-learn official documentation. *Permutation feature importance*. [Explanation and related feature constraints](https://scikit-learn.org/stable/modules/permutation_importance.html). Accessed 2026-09-09.

[^19]: scikit-learn official documentation. *Cross-validation: evaluating estimator performance*, grouped data part. [cross-validation instructions](https://scikit-learn.org/stable/modules/cross_validation.html). Accessed 2026-09-09.

[^20]: Shah RD, Peters J. *The Hardness of Conditional Independence Testing and the Generalized Covariance Measure*. [Author preprint](https://arxiv.org/abs/1804.07203), [Institutional inclusion record](https://atiro.turing.ac.uk/esploro/outputs/journalArticle/The-hardness-of-conditional-independence-testing/9915867609548). It is used to illustrate that conditional independence testing requires clear assumptions, and is not used to directly equate the prediction gain of a finite model with conditional mutual information.
