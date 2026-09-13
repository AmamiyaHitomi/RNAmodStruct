import csv
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"


def read_rows(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class HeLaModelTransferTests(unittest.TestCase):
    def test_run_status_and_population(self):
        row = read_rows(ROOT / "results" / "status" / "09c_hela_model_transfer_run_status.csv")[0]
        self.assertEqual(row["status"], "PASS")
        self.assertGreaterEqual(int(row["all_sites"]), 500)
        self.assertEqual(int(row["models_transferred"]), 5)
        self.assertLess(int(row["strict_sites"]), int(row["all_sites"]))

    def test_all_five_models_have_finite_metrics_in_both_subsets(self):
        rows = read_rows(TABLES / "09c_hela_transfer_metrics.csv")
        subsets = {row["subset"] for row in rows}
        self.assertEqual(subsets, {"all_qualifying", "strict"})
        for subset in subsets:
            models = {row["model"] for row in rows if row["subset"] == subset}
            self.assertEqual(models, {"M0", "M1", "M2", "M3", "M4"})
        for row in rows:
            for name in ["mae", "rmse", "spearman", "r2"]:
                self.assertTrue(math.isfinite(float(row[name])))

    def test_primary_comparison_is_finite_and_reported_for_both_subsets(self):
        rows = read_rows(TABLES / "09c_hela_transfer_model_comparison.csv")
        primary = [row for row in rows if row["comparison"] == "M3_vs_M1"]
        self.assertEqual(len(primary), 2)
        subsets = {row["subset"] for row in primary}
        self.assertEqual(subsets, {"all_qualifying", "strict"})
        for row in primary:
            values = [float(row[name]) for name in [
                "delta_mae_baseline_minus_augmented", "bootstrap_ci95_low", "bootstrap_ci95_high"
            ]]
            self.assertTrue(all(math.isfinite(value) for value in values))
            self.assertLessEqual(values[1], values[2])

    def test_predictions_cover_population_once_and_stay_bounded(self):
        rows = read_rows(TABLES / "09c_hela_transfer_predictions.csv")
        self.assertGreater(len(rows), 500)
        self.assertEqual(len({row["site_id"] for row in rows}), len(rows))
        for row in rows:
            for model in ["M0", "M1", "M2", "M3", "M4"]:
                self.assertGreaterEqual(float(row[f"prediction_{model}"]), 0.0)
                self.assertLessEqual(float(row[f"prediction_{model}"]), 1.0)

    def test_strict_subset_has_no_development_overlap(self):
        rows = read_rows(TABLES / "09c_hela_transfer_strict_subset_summary.csv")
        by_subset = {row["subset"]: row for row in rows}
        self.assertEqual(int(by_subset["strict"]["overlaps_development"]), 0)
        self.assertLess(int(by_subset["strict"]["sites"]), int(by_subset["all_qualifying"]["sites"]))


if __name__ == "__main__":
    unittest.main()
