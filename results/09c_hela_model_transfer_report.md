# Stage 09C HeLa frozen-model transfer report

## Scope

The five frozen HEK293T Ridge models (M0--M4) and their single frozen preprocessing encoder were applied
directly to the HeLa model population (24,960 sites with complete 201-nt sequence and valid
experimental flank reactivity, 4,687 genes). No HeLa label, feature, or
parameter was used to select or refit anything.

## Population and strict subset

| Subset | Sites | Genes | Joint groups |
|---|---:|---:|---:|
| All qualifying | 24,960 | 4,687 | 4,685 |
| Strict (no HEK293T-development gene/sequence overlap) | 22,256 | 4,048 | 4,046 |

A HeLa site is excluded from the strict subset when its analysis gene or its 201-nt sequence also appears
in the HEK293T development set used to train the frozen models (2,704
sites excluded).

## Transfer performance (MAE in proportion units; ×100 = percentage points)

| Model | Inputs | All MAE | Strict MAE | All Spearman | Strict Spearman |
|---|---|---:|---:|---:|---:|
| M0 | common background | 0.20515 | 0.20556 | 0.3622 | 0.3575 |
| M1 | background + sequence | 0.21608 | 0.21859 | 0.3231 | 0.2985 |
| M2 | background + sequence + predicted structure | 0.21891 | 0.22210 | 0.3009 | 0.2723 |
| M3 | background + sequence + experimental R_flank10 | 0.21673 | 0.21922 | 0.3231 | 0.2985 |
| M4 | background + sequence + predicted + experimental | 0.21947 | 0.22264 | 0.3009 | 0.2724 |

The primary external comparison is ΔMAE = MAE(M1) − MAE(M3): all-qualifying
-0.064 percentage points (95% CI
-0.072 to -0.056);
strict subset -0.063 percentage points (95% CI
-0.071 to -0.055).
Positive values mean lower error after adding experimental structure.

## Transfer caveats

The frozen encoder imputes HeLa's missing icSHAPE abundance field (literal `*`) with the HEK293T
development median, and one-hot categorical levels discovered only in HEK293T development are applied
as-is; HeLa rows carrying unseen `drach_subtype` motifs contribute all-zero columns for that block. Those
unseen-level counts are recorded in `09c_hela_unseen_categorical_levels.csv`. These results are direct
transfer evidence, not a HeLa-tuned model.
