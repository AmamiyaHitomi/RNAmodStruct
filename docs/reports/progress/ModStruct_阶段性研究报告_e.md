# ModStruct Phased Research Report (Phase 1 to Phase 4)

## What is this project researching?

RNA can be imagined as a long rope consisting of four letters, A, U, G, and C. It has a fixed alphabetical order and also collapses into different shapes. Certain letters are also chemically marked by cells, and these marks are RNA modifications.

ModStruct wants to answer an intuitive question:

>Do chemical labels on RNA prefer to appear in certain folding environments?

For example, if a location is more exposed and easier to be hit by chemical probes, will its m6A modification ratio be higher? Do different modifications such as m6A, m5C, m7G and Nm prefer different structural environments?

If the answer is clear, we can draw a "map of RNA modifications and structures" and even summarize the structural rules of different modifications. This was the initial vision for the project.

## Where does the data come from?

This project did not grow its own cells or conduct wet experiments, but used publicly available data from other research teams.

One type of data tells us where modifications are present in RNA, or the approximate proportion of modifications at a certain position. Another type of data uses chemical probes to measure the degree of reaction at various locations on the RNA. The degree of reaction can provide structural clues: some positions are exposed, others are affected by base pairing, protein binding, or the surrounding environment.

Structural data is not a direct picture of RNA. It's more like a record left after touching the RNA with a probe, so the measured signal isn't just affected by secondary structure. Modification data are not always in the same form: some give precise positions and proportions, some only give an enrichment interval, and some only give a list of sites. These data cannot be lumped together indiscriminately.

## Why does it take a lot of work to sort out the data in the early stage?

Different studies often use different coordinate systems. The same "100th position" may refer to the 100th position on the genome, or it may refer to the 100th position on a certain transcript; some numbers start from 0, and some start from 1. The same gene may also produce multiple different versions of RNA.

Therefore, the project first checked the data source, reference version, coordinates, positive and negative strands, transcript version and field meaning, and then aligned the modification position with the structural position. Data with insufficient structural coverage or whose locations cannot be reliably mapped will not enter the formal analysis.

This step is equivalent to confirming that the addresses on the two maps really point to the same place before merging the two maps. If the position is right or wrong, no matter how complex the statistical method is, the conclusion will be meaningless.

## How should we read the numbers in the report?The numbers in the project fall roughly into three categories: how much data was used, how big a difference was observed, and how uncertain the difference was.

| How to write the report | Popular explanation |
|---|---|
| `n=3,705 sites` | There are 3,705 modified positions that meet all criteria for this analysis. It does not mean that 3,705 independent cell experiments were performed |
| `β/SD=0.010` | When the structural signals differ by 1 standard deviation, the modification ratios differ by about 0.010 on average. The grooming scale is expressed as 0–1, so this is approximately 1 percentage point |
| `95% CI 0.007–0.013` | The data allow for an effect range of approximately 0.7–1.3 percentage points. The wider the interval, the more uncertain the estimate |
| `95% CI across 0` | The data allows for slight increases, no changes or slight decreases at the same time, and the direction of the relationship cannot be confirmed |
| `p<0.05` | Under the null effect assumption, the probability of current or more extreme data is low. It does not mean that the conclusion has a 95% probability of being true, nor does it mean that it is very effective |
| `q-value` | A correction to the p-value after checking many relationships simultaneously, used to reduce accidental positives |
| `MAE=0.20` | If you predict a modification ratio of 0–1, the prediction will differ from the measured value by about 20 percentage points on average. Smaller is better |
| `ΔMAE<0` | After adding a new set of information, the prediction error did not decrease, but increased slightly |
| `Spearman ρ=0.40` | The rankings of the two measurements are moderately consistent, but they do not give the same absolute value |

Below are the most important data for the entire project. Each line describes both the number and the question it answers.

| Research session | Specific data | What does this mean |
|---|---:|---|
| Initial data size | 4,409 associated sites for HEK293T and 25,996 for HeLa | Both cell lines have thousands to tens of thousands of analyzable sites, but truly independent experimental replicates are far fewer than the number of sites |
| Two cell lines combined | β/SD=0.01007, 95% CI 0.00713–0.01301 | A 1 SD increase in structural signal is associated with an average increase in m6A ratio of approximately 1.01 percentage points, and the relationship exists but is small |
| After subtracting known background | HEK293T correlation 0.0222;HeLa correlation 0.0141 | The residual relationship is very weak. HeLa can statistically distinguish it, HEK293T doesn't hit the same standard |
| In vivo vs. in vitro comparison | 3,705 isotopes; in vivo effect 0.01086, in vitro 0.00239 | The in vivo correlation is approximately 0.85 percentage points higher than in vitro, indicating that the cellular environment may be important |
| Change to m6A measurement technology | 2,440 HeLa sites; correlation between the two technologies ρ=0.400 | There is only moderate consistency in the ranking of the two technologies, and the structural correlation is not clear |
| Other modifications | m5C 420 sites, m7G 190 intervals, Nm 1,281 pairs; three-term q=0.749 | No clear differences were detected in the current data, but it cannot be proved that they are absolutely unrelated to the structure |
| Drug treatment | 24,898 transcripts; median control reactivity 0.23241, treatment group 0.25278; pairwise shift 0.02384 | Overall shift of structure summary signal after drug treatment; 0.02384 is structural experimental unit, not m6A changed 2.38% |

## What have we done before the first issue?

The first phase is not the starting point for the project. Before it, the project had completed a complete study using m6A: first understand the original file, then align the two sets of coordinates, screen out data that can truly be compared fairly from hundreds of thousands of candidate positions, discover relationships in HEK293T, build a model, and finally use HeLa to check whether the conclusion can be repeated in new cells. The first to fourth issues continue to ask questions based on this set of foundations.

### Step 1: Find out what each number is

`Ratio=Acov/AGcov` in the GLORI table can be understood as the estimated proportion of a site with m6A, taking the value 0–1. The project checked four tables row by row for a total of 708,103 records, and none of them did not comply with the calculation relationship. It ensures that the "1 percentage point improvement" mentioned later refers to the modification ratio, not the sequencing depth or peak fraction.

The project also confirmed that the plus-strand site is A on the genome; the minus-strand site looks like T on the genome, but still corresponds to A when read along the direction of RNA transcription. Six data files and nine reference files were checked for file size, hashing and complete readability, with 0 bad lines and structure vector length errors. Candidate data is approximately 458.7 MB. These checks indicate that the input file is complete and can be parsed correctly, but they do not mean that the original experiment is completely free of bias.### Step 2: Build a bridge between two different versions of maps

The structural data for HEK293T use the older hg19/GRCh37.74 coordinates, and the m6A data use the newer GRCh38 coordinates. The project established a conversion chain of "transcript position → old version of genome → new version of genome" and processed exon splicing and positive and negative strands.

Thirty positions covering the plus strand, minus strand, exon interiors, and splicing boundaries were selected before the official run. The center letters of the 30 RNAs are all A, and all of them get the new coordinates in the correct direction and can be converted back and forth, and the failure is 0. It demonstrates that the transformation rules pass specialized testing, but does not guarantee structural experimental coverage of every modification site.

### Step 3: Explain why there are only a few thousand sites left among hundreds of thousands

The two GLORI measurements of HEK293T had 214,715 and 215,021 unique sites, respectively, resulting in a combined total of 253,087 candidate sites.

| Screening steps | Remaining sites | What do the numbers represent |
|---|---:|---|
| GLORI merged set | 253,087 | Appears in at least one modification measurement |
| Can be mapped to icSHAPE | 20,559 | Modification and structural data can point to the same RNA location |
| GLORI detected twice | 16,161 | Two counts can be combined to estimate modification ratio |
| Meets the ±10 nt structure window | 4,409 | The structure signals on both sides are sufficient for correlation analysis |
| With full 201 nt input | 4,096 | Complete sequences and features, ready for predictive analysis |

The largest loss comes from structural coverage: 230,866 candidate sites do not have qualified icSHAPE transcript exon coverage, accounting for 91.22%; an additional 939 positions cannot switch from hg19 to hg38, and 722 positions are ambiguous. Exclusion usually means that the two experiments did not result in comparable observations at the same location and does not mean that there is no structure or no m6A at these locations.

The 4,409 associated sites come from 1,324 genes, and a single gene contributes up to 38 sites, so positions in the same gene are regarded as a group for statistics. A further 7,266 mappable positions corresponded to multiple transcript versions, up to 19, so we subsequently examined whether alternative transcript selection methods changed the conclusions.

### Step 4: First confirm whether the modification measurement is stableThe Pearson correlation of the two GLORI measurements of HEK293T was 0.9879, the Spearman correlation was 0.9728, and the median absolute difference was 0.01549. That is, the high and low rankings of the two measurements are very consistent, and the two modification ratios of a site usually differ by about 1.55 percentage points.

This only answers "can modification measurements be repeated", not "whether structure is related to modification". The latter issue needs to be modeled separately.

### Step 5: See the relationship for the first time in HEK293T

The main result is β/SD=0.01430 with a 95% confidence interval of 0.00595–0.02264. In layman's terms, a 1 standard deviation difference in structural reactivity results in an average difference in m6A proportions of approximately 1.43 percentage points; the data allowable range is approximately 0.60–2.26 percentage points. Seven prespecified analytical approaches all yielded positive results, ranging from approximately 1.02–2.05 percentage points.

This relationship can be measured, but the amplitude is not large, and it is an observational comparison between natural sites. It cannot be said that "artificially changing the structure will definitely increase m6A by 1.43%." The model R²=0.1742 means that the structure, sequence region and all other variables jointly explain about 17.4% of the site differences. It cannot be said that the structure alone explains 17.4%.

### Step 6: Check whether the structure can help predictions

The project used 4,096 sites to train the model, 3,271 for development, and 825 reserved for internal testing. When dividing, ensure that there are no identical genes, the same 201 nt sequence, or repeated groups connected by the two on both sides to prevent the model from recognizing the trained close relative samples during testing.

Using only the M1 model of background and sequence, the test MAE=0.20546, that is, the average difference between the predicted and measured modification ratios is about 20.55 percentage points. After adding the computer prediction structure, M2 MAE=0.21193, the error increased by about 0.65 percentage points. Adding the M4 MAE=0.21188 of the experimental structure, the improvement is only 0.0052 percentage points relative to M2, and the confidence interval spans 0. It is impossible to confirm that this difference is not a random fluctuation. Development efforts focused on attempting a more complete experimental structure profile also did not outperform models using only background and sequence.

So "structure and modification change together" and "structure makes model predictions more accurate" are two issues. There can be a relationship between structure and modification, but the information carried by the structure may also have been included in sequence and region variables.

### Step 7: Use HeLa to do the real external examThe HEK293T internal test set was later repeatedly viewed in preprocessing revisions and supplementary analyses, and the project therefore froze the rules and models for direct use on HeLa data that had never been involved in development.

HeLa started with 164,761 GLORI merged sites and ended up with 25,996 association sites and 24,960 model sites, with an association set from 4,829 genes. The main association was 0.00947/SD with a 95% confidence interval of 0.00633–0.01261, i.e., when structures differed by 1 SD, m6A proportions differed on average by about 0.95 percentage points. It's smaller than HEK293T's ~1.43 percentage points, but in the same direction; all seven preset checks also remain positive. This suggests that the "small positive association" was repeated in another cell line.

The frozen model was applied directly to 24,960 HeLa sites; after excluding all sites that are isogenic or iso-sequential to the HEK293T development set, the strict subset remained at 22,256. In the entire HeLa, M1 MAE=0.21608, and the M3 MAE=0.21673 added to the experimental structure, the error worsens by about 0.064 percentage points; the strict subset also deteriorates by about 0.063 percentage points. HeLa was not used for model selection or parameter retuning, so this is a more rigorous external test.

### Full conclusion before the start of the first issue

The preliminary research has completed the entire link of "understanding the data - verifying coordinates - quantitative screening - HEK293T discovery - leak prevention modeling - HeLa external review", and saved frozen configurations, result hashes and automated checks. The HEK293T correlation analysis failed none of the 14 tests and the overview chart passed all 21 quality checks. This shows that the computational process can be stably reproduced, but does not mean that biological causation has been proven.

> There is a repeatable but small positive relationship between local experimental structure and m6A ratio, which is about 1 SD difference in structure and 1 percentage point difference in modification ratio; however, after sequence and background information are entered into the model, experimental structure does not show stable additional predictive value.

The next four periods will focus on this sentence to check whether the answer still holds true after changing statistical methods, structural conditions, measurement techniques, modification types and disturbance backgrounds.

## Issue 1: Conduct a comprehensive physical examination of this conclusion

The purpose of the first period was to check whether the initial results were caused by just one analysis method.

We have summarized the results for HEK293T and HeLa. The combined effect is approximately 0.010,That is, for every one standard deviation increase in structural indicators, the m6A proportion changes by about one percentage point on average. This relationship is not significant, but the directions of the two cell lines are basically the same.

More specifically, the fixed effects pool is 0.01007 and the 95% confidence interval is 0.00713–0.01301. After scaling, the average relationship falls approximately between 0.71–1.30 percentage points. The comparison here is between naturally occurring different sites, so it cannot be interpreted as "if the structure is artificially increased by 1 SD, m6A will increase by 1%."

Subsequently, various analysis methods were changed: examining associations after subtracting sequence and background information, examining the degree of uncertainty in predicted structures, changing the selection method of transcript versions, expanding the observation window from 201 nt to 401 nt, and using models that can handle complex nonlinear relationships.

The overall answer from these checks does not change: small correlations are still seen in several analyses, but the experimental structure does not steadily improve external predictive power. More complex models predict better overall, without demonstrating additional help from the structural information itself.

Taking HeLa external data as an example, in the 401 nt window analysis, the MAE of the model using only background and sequence was 0.22758, and after adding both predicted and experimental structures, the MAE was 0.22875. The error increases by about 0.00117, which is equivalent to about 0.12 percentage points. The complete model MAE of the nonlinear model is 0.18942, which is lower than Ridge's 0.21463; however, compared within the nonlinear model, the error after adding the experimental structure still increases by about 0.00169. This shows that the nonlinear algorithm itself is more suitable for some data relationships, but it does not prove that the experimental structure provides new information.

The first issue tells us that the initial result does not disappear completely simply by switching to a common method, but it is also not a strong and easily predicted relationship.

## Issue 2: Change the environment and change the measurement method and watch it again

The first phase mainly involves changing the analysis methods on the same batch of data. The second issue further asked: If the experimental environment is changed or the measurement technology is modified, will the answer remain the same?

The project began by comparing in vivo and in vitro structural signals from the same m6A sites in HEK293T. There ended up being 3,705 locations for a fair comparison. The correlation was approximately 0.0109 in vivo and 0.0024 in vitro; the confidence interval for the difference did not cross zero.

The difference between the in vivo and in vitro effects was specifically 0.00848, with a 95% confidence interval of 0.00152–0.01519. Switch to proportional language,That is, the in vivo correlation is on average about 0.85 percentage points higher than in vitro , and the data allow for a difference of about 0.15–1.52 percentage points.

This suggests that the relationship between m6A and structure may depend on the intracellular environment. RNA-binding proteins, ions, molecular crowding and various regulatory processes exist in cells, and the in vitro environment does not fully retain these factors. However, this comparison cannot yet tell us which factor is at play.

The project also compared the original GLORI data with another m6A measurement technique, SAC-seq. A total of 2,440 HeLa sites were able to enter the analysis simultaneously. The values ​​measured by the two techniques were moderately correlated, but the structural correlations were small, with confidence intervals spanning zero; structure did not improve the predictive performance of SAC-seq results.

Specifically, the Spearman correlation for the two modification measures was 0.400. The effect of experimental structure was 0.00407/SD with a 95% confidence interval of -0.00436 to 0.01249 for GLORI and 0.00121/SD with a 95% confidence interval of -0.00465 to 0.00707 for SAC-seq. Because both intervals cross zero, this batch of matched data does not confirm a stable positive or negative structural effect.

Therefore, the most important findings of Issue 2 are:

> The relationship between modification and structure will be affected by experimental conditions and measurement techniques, and the results in a certain data set cannot be directly regarded as a universal rule.

The Human RNA MaP and PARS data originally considered were not forcibly used because no quantitative m6A data matching the conditions could be found. This is not a task failure, but the result of enforcing a data quality threshold.

## Issue 3: Do other modifications have their own structural preferences?

The third phase returned to the project's original vision of multi-modified maps. The project analyzed m5C, m7G and Nm.

The three types of data are in different forms. The m5C data is closer to the modification ratio, the m7G data is mainly the enrichment degree in an interval, and the Nm is closer to the modification site list. They are like three maps with different accuracy, and the numbers cannot be added together directly. So the project analyzed each modification separately and then compared the overall evidence.

The result is as follows:

- m5C has 420 qualified sites in HeLa with an effect of -0.00349/SD and a 95% confidence interval of -0.01990 to 0.01292, with no clear structural association detected.
- m7G has 190 qualified intervals in HeLa,The interval enrichment effect was -0.02636/SD, the 95% confidence interval was -0.18790 to 0.13519, and no clear interval-level association was detected.
- Nm formed 1,281 sets of modification sites and similar controls in HeLa, with an average structural difference of -0.00451 and a 95% confidence interval of -0.02014 to 0.01217, and the structural difference between the two was also unclear.

The confidence intervals for the three main tests all included zero, and none passed the multiple testing criterion. This does not prove that these modifications are absolutely independent of structure, as sample size, measurement error, and feature selection all affect detection power. It can only show that using current public data and current designs, no clear and reliable differences have been found.

Therefore, the third phase has completed the reconnaissance, review, alignment and first round of comparison of multi-modification data, but it has not been able to draw a credible map of "what structures different modifications prefer".

## Issue 4: After changing the system, will the structure change accordingly?

The first three parts focus on observing co-variations in natural data. The fourth issue uses a set of drug treatment data to take a step closer to the directionality issue.

STM2457 inhibits METTL3 in association with m6A addition. Published data compare STM2457-treated and control groups in the same HEK293T background. The project paired 24,898 transcripts and found a systematic shift in average structural reactivity after treatment, with a median change of approximately 0.0238.

The median average transcript reactivity was 0.23241 for the control group and 0.25278 for the treatment group. The algorithm for pairwise shift 0.02384 is to first calculate "treatment minus control" for each homonymous transcript and then take the median of all differences, so it does not need to be equal to the direct subtraction of the two group medians. It uses normalized SHAPE reactivity units, with higher values ​​generally indicating greater access of the probe to those RNA locations, but the signal can also be affected by factors such as protein binding.

This suggests that the RNA structure summary signal does change after drug treatment occurs. But it’s not possible to say yet that “m6A reduction leads to structural changes.” Reasons include:

- What the project obtains is the summary result of each transcript, not the changes of each m6A site in the same sample.
- Drugs may affect structure through RNA content, protein binding, or other pathways.
- Lack of controls such as rescue, inactive drugs, and independent data review.
- 24,898 transcripts represent many subjects and do not equate to 24,898 independent experiments.So what arrives in the fourth period is the "directional clue": there is inhibitory processing first, and then structural changes are seen. It has not yet reached the point of "proving cause and effect", let alone "explaining the specific mechanism".

The project did not train Transformer or perform large-scale parameter adjustment in this phase because the existing data lacked sufficient modification types, processing conditions, and independent test sets. At this time, increasing model complexity cannot fill the gap in experimental design.

## What do the four issues combined mean?

The entire project can be summarized in one story line:

```text
initial guess
Different RNA modifications may prefer different structural environments
        ↓
m6A baseline study
Small structural correlations found, but structure did not steadily increase predictive power
        ↓
Methods, environment and technical review
The correlation will change with conditions and is not a strong, fixed rule.
        ↓
Explore multiple modifications
No clear modification-specific structural fingerprint has yet been obtained
        ↓
drug disturbance
The structure shifted after treatment, but the reason has not yet been proven
```

The findings were not exactly what was initially expected. We don't get a brightly colored map with clear structural preferences for each modification. The answer was more cautious:

> There may be a weak, context-dependent relationship between RNA modification and structure; how much of this is information about the structure itself and how much comes from sequence, measurement methods, and cellular context remains to be determined.

This is still a valid scientific research result. More complex models that have been tried so far do not address structural gain; this makes comparable data and independent validation a higher priority, but cannot rule out the possibility that other structural features or models may provide gains in the future.

## Looking back, what is the idea and logic of the entire project?

The analysis process of the current project is relatively solid, but the main issues need to be addressed more focused. m6A is the most well-studied for association and prediction, with other modifications and drug treatments providing extended clues and yet to collectively demonstrate a complete set of mechanisms. This is an evaluation of existing results, without adding new experiments or re-running all analyses.

The best thing about the project is to first check whether the data corresponds correctly, then do the analysis, and finally use another cell background to check the results. It also separates two questions: "Is there a relationship between structure and modification ratio?" and "Can predictions be more accurate after knowing the structure?" There is a slight relationship, but there is no stable improvement in prediction accuracy, which can happen at the same time. The subsequent changes of windows, models, and conditions are to check how well this answer can withstand the test.

What needs to be made most clear is that the original question is different from the question being answered now.Initially, I wanted to know "where modifications are more likely to occur, and whether different modifications prefer different structures." The most complete analysis at present is "where modifications have been detected and structures can be measured, what is the relationship between the modification ratio and the structure." Knowing how the proportion of labels changes in places that already have labels cannot explain why there are no labels elsewhere. To answer the latter question, credible unmodified or low-modification controls are required.

Not every completion of the four phases of work will bring you closer to proving cause and effect:

- The early stage and the first period are the main subjects to study whether the correlation is reliable and whether the structure can improve predictions.
- The second phase checks whether the results can be seen again after changing the environment and technology.
- The third phase explores other modifications, but the ratio, enrichment interval and site list cannot be directly ranked into a structural preference list.
- The fourth phase observes structural changes after medication, but it has not been proven that the changes are caused by the reduction of m6A.

Therefore, completing four phases means that this round of planning is completed, but it does not mean that all the laws and mechanisms originally envisioned have been found.

Another boundary comes from sample screening. HEK293T initially had approximately 253,000 merged sites, with 4,409 eventually entering the association analysis. Screening is to ensure data quality, but the answer first applies to these positions where "modifications can be detected and structures can be measured" and cannot automatically represent all RNA positions.

Results after changing measurement techniques should also be interpreted with caution. Comparisons between in vivo and in vitro directly tested the differences and provided clues that environment may be important; matching analysis of GLORI and SAC-seq did not confirm a clear association, indicating that support for cross-technical replication was insufficient. Just because one result is significant and the other is not, it cannot be said directly that the effects obtained by the two techniques are indeed different.

The best questions to focus on right now are:

> How reliable is the relationship between local structure and m6A ratio? Now that we know the sequence and background, can the structure help us predict more accurately on new data?

The current answer is: a small correlation can be seen, but no stable additional predictive help within the range of data, structural indicators and models examined. This does not mean that structures have no biological role, nor does it mean that all structural information is useless.

The next step should be to complete this main line first, with m6A as the center, and put multi-modification and drug results in an expanded position. If you continue to look for local patterns in the context of certain regions or sequences, you should specify a small number of questions in advance and test them with new independent data; if you cannot see them repeatedly, keep the results as they are. What is most worthy of improvement at present is the consistency between "what is asked, what is actually measured, and what can be said in the end".

## What is the most worthwhile next step?The existing four issues can be compiled into a complete research version. The main body should be organized around m6A: reporting small associations first, then presenting various robustness tests, in vivo and in vitro differences, cross-technology comparisons, and prediction boundaries. Other modifications and pharmacological treatments are suitable as extended results to illustrate where the current evidence goes and where it stops.

The questions most worthy of answering in the new round of research are:

> Under what conditions does this small correlation exist, and can it be seen again in a new batch of data that was never used in the development analysis?

A small number of meaningful classifications can be specified in advance, such as 5′UTR, CDS, and 3′UTR of RNA, different sequence motifs, and different coverage levels. The complete structure curves around the modification site are then compared and a check is made to see whether a certain type of local pattern recurs in independent data.

If you continue to do drug perturbation, you should prioritize obtaining data for each biological replicate and try to find site-by-site m6A changes, structural changes, and more complete controls in the same sample. This way it can be determined whether the structural changes are actually caused by changes in m6A.

There is no need to rush into adding more polish, building a large website, or training a large model. What the project currently lacks most is comparability and independent verification between different data. Only after these problems are solved can complex models and structural maps have a reliable foundation.

## Conclusion of the current stage

ModStruct has completed a complete round of research from data collection, m6A mainline, method physical examination, cross-condition and cross-technology review, to multi-modification exploration and drug perturbation analysis.

What can be said with confidence now is that there is a small, condition-sensitive, observational relationship between m6A ratio and local experimental structural signal. It is not yet possible to say that structure consistently improves modification predictions, that various RNA modifications have shown clear distinct structural preferences, or that m6A reduction has been shown to lead to changes in RNA structure.

This stage can be concluded and a full report produced. The next stage should focus on finding local structural regularities that can be repeated in independent data, while using more appropriate experimental units and control designs, so that the project can gradually move from "seeing relationships" to "explaining relationships."

## Data location

The original concept, four-phase plan and results of each phase of the project are recorded in the following documents:

- [Original project idea](../../research_notes/chat_idea_organized.md)
- [Fourth Phase Restart Scope List](../../research_notes/ModStruct_Restart Scope List.md)- [GLORI field semantic verification](../../../metadata/provenance/01_glori_field_semantics.md)
- [Data source and file integrity verification](../../../metadata/provenance/source_verification.md)
- [HEK293T coordinates and initial intersection report](../../../results/reports/05_hek293t_initial_intersection_report.md)
- [HEK293T analysis data report](../../../results/reports/06_hek293t_analysis_dataset_report.md)
- [HEK293T Association Validation Report](../../../results/reports/07_hek293t_association_validation_report.md)
- [HEK293T prediction validation report](../../../results/reports/08_hek293t_prediction_validation_report.md)
- [HeLa Analysis Data Report](../../../results/reports/09_hela_analysis_dataset_report.md)
- [HeLa Association Review Report](../../../results/reports/09b_hela_association_report.md)
- [HeLa frozen model transfer report](../../../results/reports/09c_hela_model_transfer_report.md)
- [Phase 1 Method Expansion Summary](../../../results/reports/phase1_summary_report.md)
- [Phase 2 Completion Report](../../../results/reports/22_phase2_completion_report.md)- [Phase 3 Completion Report](../../../results/reports/27_phase3_completion_report.md)
- [Phase 4 Completion Report](../../../results/reports/31_phase4_completion_report.md)