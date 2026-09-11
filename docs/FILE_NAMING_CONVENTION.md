# File naming convention

All pipeline-owned files use the project-stage number from the research plan as their first token.

## Stage prefixes

| Prefix | Project stage | Content |
|---:|---|---|
| `01_` | GLORI semantics | Field definitions and evidence |
| `02_` | Reference preparation | Reference downloads, manifests, and compatibility checks |
| `03_` | Input audit | Integrity, schema, missingness, and range audits |
| `04_` | Coordinate validation | Small stratified mapping checks and assertions |
| `05_` | Preliminary intersection | Full mapping status, joined tables, attrition, and summaries |
| `06_` | Analysis dataset build | Frozen configuration, isoform resolution, final tables, and QC |
| `07_` | Association analysis | Descriptive QC, clustered inference, sensitivity checks, and figures |
| `08_` | Prediction modeling | Frozen group split, Ridge ablation, held-out evaluation, and model artifacts |
| `09_` | HeLa replication | HeLa dataset build, predicted structure, association re-estimate, and frozen-model transfer |

## Directory roles

- `src/`: executable code only. Names start with the stage and an action verb.
- `metadata/audit/`: machine-readable stage-03 audit evidence only.
- `metadata/coordinate/`: machine-readable stage-04 coordinate validation evidence only.
- `data/interim/`: derived row-level tables; every filename carries its producing stage.
- `data/final/`: analysis-ready datasets; every filename carries its producing stage and population.
- `results/`: compact result summaries and human-readable reports; filenames carry their producing stage.
- `results/logs/`: one log per executable stage, with the same stage/action stem as its script.
- `tests/`: tests named after the stage or shared module they validate.

Shared library code that is not an executable stage may omit a prefix, for example
`src/pipeline_common.py`. Immutable downloaded source filenames remain unchanged so they can be compared
directly with their official accession records.

## Pattern

Use lowercase snake case:

```text
<stage>_<action-or-scope>_<data-object>.<extension>
```

Do not use vague names such as `final.csv`, `output.csv`, `results2.csv`, or reuse the same unqualified
basename in different stages.
