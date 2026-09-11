# Stage 07 figure contract

Core conclusion: quantify the direction and magnitude of the prespecified association between local
icSHAPE reactivity and GLORI m6A proportion, and show whether that estimate persists under predefined
robustness checks.

Results-level question: within repeatedly detected HEK293T GLORI sites with adequate experimental
structure coverage, how is local icSHAPE reactivity associated with the pooled m6A proportion?

Figure archetype: quantitative grid.

Target/output: report-ready, double-column scientific figure; editable SVG is primary, with PDF and TIFF
companions.

Backend: Python/matplotlib exclusively.

Final size: 180 mm × 105 mm; all rendered text at least 5 pt. Panel a spans both rows as the hero
description; panels b and c are stacked at right.

Panel map:

- a — raw population density and binned trend; primary descriptive evidence.
- b — adjusted primary coefficient and predefined sensitivity estimates; robustness/boundary evidence.
- c — ±50 nt mean reactivity profiles by frozen m6A-ratio quartile with pointwise cluster-bootstrap
  intervals; spatial descriptive context rather than an additional primary test.

Evidence hierarchy: panel a is the hero description, panel b carries inferential robustness, and panel c
shows where measured signal is available around the site.

Statistics: OLS coefficient, gene-cluster robust CI, gene-cluster bootstrap CI, exact sample/gene counts,
and BH adjustment across the seven sensitivity analyses. Profile bands are pointwise rather than
simultaneous confidence bands.

Source data: stage-07 machine-readable association, sensitivity, binned-trend, and positional-profile
tables.

Image integrity: vector-native statistical plots; no raster image manipulation. Dense raw observations
may be represented by hexagonal binning without subsampling.

Reviewer risks: selected-site/collider bias, residual nonlinearity, bounded outcome with heteroscedastic
errors, cross-study batch confounding, gene-level dependence, missing-structure selection, isoform
ambiguity, and inability to infer direction or causation.
