"""Re-run stages 12--17 and require byte-identical machine-readable outputs."""
from __future__ import annotations

import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from phase1_common import sha256_file, write_rows

ROOT = Path(__file__).resolve().parents[1]
COMMANDS = [
    "12_run_meta_analysis.py", "13_run_conditional_dependence.py", "14_run_structure_entropy.py",
    "15_run_isoform_sensitivity.py", "16_run_window_401_sensitivity.py", "17_run_nonlinear_models.py",
]
OUTPUTS = [
    "results/tables/12_meta_analysis.csv",
    "results/tables/13_conditional_dependence.csv", "results/tables/13_cross_fitted_residuals.csv.gz",
    "results/tables/14_structure_entropy_associations.csv",
    "results/tables/15_isoform_sensitivity_associations.csv", "results/tables/15_isoform_site_estimates.csv.gz",
    "data/interim/16_401nt_vienna_cache.csv", "data/final/16_hek293t_401nt_dataset.csv.gz",
    "data/final/16_hela_401nt_dataset.csv.gz", "results/tables/16_window_401_attrition.csv",
    "results/tables/16_window_401_development_cv.csv", "results/tables/16_window_401_external_metrics.csv",
    "results/tables/17_nonlinear_development_cv.csv", "results/tables/17_nonlinear_external_metrics.csv",
    "results/tables/17_grouped_permutation_importance.csv", "results/tables/17_retraining_ablation.csv",
]
TABLE = ROOT / "results/tables/phase1_reproducibility_check.csv"
STATUS = ROOT / "results/status/phase1_reproducibility_status.csv"


def main() -> None:
    started = datetime.now(timezone.utc)
    paths = [ROOT / value for value in OUTPUTS]
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing original outputs: {missing}")
    before = {path: sha256_file(path) for path in paths}
    timings = {}
    for script in COMMANDS:
        began = time.monotonic()
        subprocess.run([sys.executable, str(ROOT / "src" / script)], cwd=ROOT, check=True)
        timings[script] = time.monotonic() - began
    rows = []
    for path in paths:
        after = sha256_file(path)
        rows.append({"path": path.relative_to(ROOT).as_posix(), "sha256_before": before[path],
                     "sha256_after": after, "verdict": "MATCH" if before[path] == after else "MISMATCH"})
    write_rows(TABLE, rows)
    passed = all(row["verdict"] == "MATCH" for row in rows)
    write_rows(STATUS, [{"status": "PASS" if passed else "FAIL", "method": "deterministic_full_rerun",
                         "verdict": "REPRODUCIBLE" if passed else "NOT_REPRODUCIBLE",
                         "started_utc": started.isoformat(), "finished_utc": datetime.now(timezone.utc).isoformat(),
                         "artifacts_compared": len(rows), "commands": len(COMMANDS),
                         "duration_seconds": sum(timings.values()),
                         "command_timings_seconds": ";".join(f"{key}={value:.3f}" for key, value in timings.items())}])
    if not passed:
        raise AssertionError("One or more phase-1 machine outputs changed on deterministic re-run")


if __name__ == "__main__":
    main()
