# ModStruct Popular interpretation: What exactly is this project doing?

> This article is the "vernacular version" of the "ModStruct Undergraduate Research Plan", using life metaphors to explain the entire project. If you want to know the rigorous details, please read the original plan; if you want to understand "what this job is for" first, this article is enough.

---

## 1. One sentence version

**We need to figure out: How strong is the relationship between the "shape" of the RNA molecule and the m6A chemical label on it - especially how much of the relationship remains after excluding other interfering factors. **

---

## 2. Add a little background first (don’t be afraid, it’s very short)

**What is RNA? **
There is a type of molecule in cells called RNA, which you can think of as a "letter chain" with only 4 letters (A, U, G, C). This chain does not lie straight, it folds on itself to form a three-dimensional shape - just like a rope winding itself into a ball, where it is rolled and where it is unrolled. This shape (technically called "structure") affects the efficiency of the RNA's work.

**m6A What is it? **
Scientists have discovered that certain "letter A" positions on the RNA chain will be labeled by cells with a small chemical label (methyl group), which is m6A. This label is very important. It affects whether RNA is stable or destroyed, and whether it is translated into protein. It is one of the hottest research directions now.

**Key Background:**
Scientists have long known that m6A and RNA structures are "related" to each other - the shape affects whether the label sticks on, and the label will in turn change its shape. Therefore, "the two are related" is not something new to be discovered by this project. This must be made clear first. Don't treat old conclusions as new discoveries.

---

## 3. So what exactly do you want to do with this project?

For example:

> Everyone knows that "the more expensive the house, the better the decoration is generally." However, housing prices are also affected by location, area, and age of the property.
> The really interesting question is: **After knowing the location, area, and age of the house, can "decoration style" additionally explain the difference in housing prices? **

Translated into this project:

- "House price" = m6A modified amount
- "area/lot/age" = RNA sequence itself (alphabetical order), gene region, measured quality
- "Decoration style" = RNA structure measured experimentally

The project will answer three specific questions:

1. **Are they related? ** For those sites that have been confirmed to contain m6A, is there any correlation between the "how much" modification and the RNA structural signal around it?2. **Does the experimental structure have additional value? ** On the basis of "only using sequences" and adding "experimentally measured structural data", will the ability to predict modification levels be improved?
3. **Is cell replacement still valid? **Do the patterns found in HEK293T cells hold in HeLa cells (another cell line)?

The second question is the most valuable. It is equivalent to asking: Is the information of "experimental structure" something **extra** useful outside of the sequence?

---

## 4. Where does the data come from? Want to do your own experiment?

**There is no need to do experiments at all**, all data published on the Internet are used. Two types of data come from different public databases (GEO):

| Data | Function | Source cells |
|---|---|---|
| GLORI data | Measuring "how much" m6A modification | HEK293T (main) + HeLa (check) |
| icSHAPE data | Measure RNA structure (measured experimentally) | HEK293T (main) + HeLa (review) |

The core engineering work to be done in the project is to align these two types of data at the "same position", put them together into a large table, and then perform analysis on this table.

**This step is also the most dangerous step** (see "Difficulties" below).

---

## 5. How to do it? (Method, speak human language)

### 5.1 "Audit" the data first (the first two weeks are the most important)

Instead of just downloading it and using it, you need to first find out what is written in each file: What does one line mean? Which column is the modification ratio? How are the coordinates calculated? This step is called "auditing". It is repeatedly emphasized in the plan: **No analysis is allowed before you understand the meaning of the data**.

### 5.2 Do correlation analysis (statistics)

Use a simple statistical model to see if there is any relationship between "structural signal" and "modification amount", how strong the relationship is, and how big the error is. At the same time, it is necessary to "control" a bunch of interfering factors (sequence characteristics, gene regions, measurement coverage, etc.) to see if the relationship is still there after eliminating these.

### 5.3 Use machine learning to do "ablation experiments" (distinguish credit)

This is the neatest part. Build 5 models (M0 to M4) and add things to them step by step:

| Model | What was added | Want to answer |
|---|---|---|
| M0 | Background information only | How much can the area and measurement itself explain? |
| M1 | + Sequence | How strong is the sequence itself? |
| M2 | + Calculate the predicted structure | Is the predicted structure useful?|
| M3 | + Experimentally tested structure | Does the experimental structure have additional value? |
| M4 | All required | Based on the predicted structure, is there any gain in the experimental structure? |

**The core comparison is M4 vs. M2**: If M4 predicts more accurately, it means that the experimental structure does bring additional information. This step can separate and calculate the "credit of the sequence" and the "credit of the experimental structure".

The simplest model (Ridge regression, gradient boosting tree) is used, **no deep learning, no GPU, no neural network training** - because it is not necessary, the simple model is enough to answer the question "Where does the information come from?"

### 5.4 Review of changing cell lines

The rules obtained on HEK293T were transferred to HeLa data intact for verification. If it still holds true in the past, the conclusion will be more convincing; if it doesn't hold in the past, write it down truthfully.

---

## 6. How to arrange time (12 weeks)

| Stage | What to do |
|---|---|
| Week 1–2 | Download data, audit files, and see if they can be aligned (the most critical) |
| Week 3–4 | Build joint data tables, do quality control, and **freeze the plan** (no changes are allowed later) |
| Weeks 5–7 | Correlation Analysis + Machine Learning Models |
| Week 8 | Come up with the "minimum usable version" (single cell line, reproducible) |
| Weeks 9–10 | Review with HeLa data |
| Week 11–12 | Compile reports and defense materials |

---

## 7. Where is the difficulty? (This is where it really gets stuck)

**Conclusion first: the difficulty is almost not in "writing code", but in "judging whether the data is correct" and "controlling yourself not to cheat". **

1. **Data alignment is most likely to be completely wrong. ** The two types of data come from different people, different years, and different versions. If there is only one mistake in the coordinate rules, gene versions, and positive and negative chains, all subsequent analyzes will be completely wrong - and the code will be run and the pictures will be shown, and the mistake will be silent.

2. **The true status of the data is not yet known. ** When the plan was written, the file was not checked field by field. For example, "whether this column modifies the proportion or something else" is not yet sure. This can only be done by actually downloading the file, opening it, and checking it column by column. You cannot rely on guessing.

3. **Beware of "cheating just to look good". ** For example, repeatedly adjust parameters until p<0.05, select the best random seeds, allow the same gene to be used for both training and testing, and only display good-looking sites. None of these ** are due to "technical incompetence",It's all about "do you want to restrain yourself"**.

4. **Keep the boundaries of interpretation. ** Some things should not be said nonsense: structural signal ≠ pairing probability, "controlled covariates" ≠ "eliminated all interference", the gain of predicted structure can only be called "representation gain" and cannot be touted as a new discovery.

5. **AI cannot help with these tasks, and is even more dangerous. **AI can help you write code, check information, draw pictures, and write reports, but it will not actively stop. You tell it to "tweak it to make it look better" and it p-hacks you all the way until it produces a beautiful but wrong result.

---

## 8. What should be handed over in the end?

- a research report;
- 5 main pictures (data source, correlation curve, effect forest diagram, model comparison, review results);
- Complete analysis code + data dictionary + filtered records;
- A "limitation list" (an honest statement of the limits within which the conclusion holds true).

---

## 9. What is the value of this project (non-bragging version)

It doesn't pursue cool methods, nor does it promise "first discovery" or "guaranteed publication." Its value lies in:

> **Use public data + rigorous statistics + simple machine learning to reproducibly answer a specific question - "Can the experimentally measured RNA structure help us predict the m6A modification level in addition to the sequence?" **

And the plan is very clear: **No matter whether the result is "with gain", "no gain" or "only with gain in a small area", it is considered a qualified completion**. The quality of the project depends on whether the evidence is reliable, not on whether the conclusion is good or not.

---

## 10. One sentence summary

**This is a well-designed undergraduate bioinformatics topic: no experiments or large models are conducted, and the "additional predictive value of RNA experimental structure for m6A modification" is quantified by aligning two types of public data + rigorous statistics + ablation comparison. The difficulty is not in programming, but in data auditing, judgment and maintaining the discipline of "not cheating". **