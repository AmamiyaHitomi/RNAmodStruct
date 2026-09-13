# ModStruct first phase operation instructions

The first phase is an exploratory expansion after freezing the HEK293T–HeLa main line, and does not redefine the analysis conclusions of the 06–10 phases.

## Running order

| Sequence | Command | Content |
|---:|---|---|
| 1 | `python src/12_run_meta_analysis.py` | Two-cell line fixed/random effects meta-analysis |
| 2 | `python src/13_run_conditional_dependence.py` | Cross-fitting residual-on-residual conditional dependence test |
| 3 | `python src/14_run_structure_entropy.py` | Structural entropy correlation and BH correction |
| 4 | `python src/15_run_isoform_sensitivity.py` | Four isomer mapping sensitivity calibers |
| 5 | `python src/16_run_window_401_sensitivity.py` | Full 401 nt sequence, icSHAPE and ViennaRNA sensitivity |
| 6 | `python src/17_run_nonlinear_models.py` | Ridge/HGB, nested ablation and grouped permutation importance |
| 7 | `python src/phase1_verify_reproducibility.py` | Deterministic full-process rerun and byte-by-byte comparison of machine products |
| 8 | `python src/phase1_finalize.py` | Phase 1 summary, statistical fallacy scan and independent checklist |

Execute using the project environment, for example:

```text
E:\ancd\envs\my_pytorch\python.exe src/12_run_meta_analysis.py
E:\ancd\envs\my_pytorch\python.exe -m pytest tests -q
```

## Interrupt recoveryPhase 16: Append writes every 50 new unique sequences.
`data/interim/16_401nt_vienna_cache.csv`. After interruption, restart with the same command to reuse the existing hash, and the completed sequence will not be refolded.

## Products and frozen boundaries

- Summary report: `results/reports/phase1_summary_report.md`
- Reproducibility check: `results/tables/phase1_reproducibility_check.csv`
- Statistical fallacy scan: `results/tables/phase1_statistical_validation.csv`
- Phase 1 independent manifest: `metadata/manifests/release/phase1_release_manifest.csv`
- Phase 1 completion status: `results/status/phase1_completion_status.csv`

Issue 1 listings are for extended product use only; the content and hash of the original 06–10 file must not be changed. A new file entry can be added to the repository-level release manifest, but the path, size, and SHA-256 of the existing baseline entry must remain unchanged.