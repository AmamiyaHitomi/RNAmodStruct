import csv
import math
import unittest
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"


def read_rows(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class PredictionModelingTests(unittest.TestCase):
    def test_run_status_and_population(self):
        row = read_rows(ROOT / "results" / "status" / "08_hek293t_prediction_modeling_run_status.csv")[0]
        self.assertEqual(row["status"], "PASS")
        self.assertEqual(int(row["population_sites"]), 4096)
        self.assertEqual(int(row["models_fit"]), 5)
        self.assertEqual(int(row["models_pending"]), 0)
        self.assertEqual(int(row["bootstrap_successful"]), 2000)
        self.assertEqual(row["test_status"], "REUSED_AFTER_PREPROCESSING_CORRECTION")

    def test_split_has_no_group_gene_or_sequence_leakage(self):
        rows = read_rows(TABLES / "08_hek293t_split_audit.csv")
        overlap = next(row for row in rows if row["split"] == "overlap_audit")
        self.assertEqual(int(overlap["gene_overlap"]), 0)
        self.assertEqual(int(overlap["sequence_overlap"]), 0)
        self.assertEqual(int(overlap["joint_group_overlap"]), 0)

    def test_model_matrix_is_explicit(self):
        rows = read_rows(TABLES / "08_hek293t_test_metrics.csv")
        status = {row["model"]: row["status"] for row in rows}
        self.assertEqual(status, {model: "FIT" for model in ["M0", "M1", "M2", "M3", "M4"]})

    def test_selected_alphas_are_inside_search_grid(self):
        config = yaml.safe_load((ROOT / "config" / "08_prediction_modeling.yaml").read_text(encoding="utf-8"))
        grid = [float(value) for value in config["learner"]["alpha_candidates"]]
        rows = read_rows(TABLES / "08_hek293t_test_metrics.csv")
        self.assertTrue(all(min(grid) < float(row["alpha"]) < max(grid) for row in rows))

    def test_primary_comparison_is_finite_and_consistent(self):
        rows = read_rows(TABLES / "08_hek293t_primary_model_comparison.csv")
        self.assertEqual(len(rows), 4)
        row = next(item for item in rows if item["comparison"] == "M4_vs_M2")
        values = [float(row[name]) for name in [
            "delta_mae_baseline_minus_augmented", "bootstrap_ci95_low", "bootstrap_ci95_high"
        ]]
        self.assertTrue(all(math.isfinite(value) for value in values))
        self.assertLessEqual(values[1], values[0])
        self.assertGreaterEqual(values[2], values[0])

    def test_predictions_cover_test_once_and_stay_bounded(self):
        rows = read_rows(TABLES / "08_hek293t_test_predictions.csv")
        self.assertGreater(len(rows), 700)
        self.assertEqual(len({row["site_id"] for row in rows}), len(rows))
        for row in rows:
            for model in ["M0", "M1", "M2", "M3", "M4"]:
                self.assertGreaterEqual(float(row[f"prediction_{model}"]), 0.0)
                self.assertLessEqual(float(row[f"prediction_{model}"]), 1.0)

    def test_all_sites_have_a_persisted_group_assignment(self):
        rows = read_rows(TABLES / "08_hek293t_group_assignments.csv")
        self.assertEqual(len(rows), 4096)
        self.assertEqual(len({row["site_id"] for row in rows}), 4096)
        self.assertEqual({row["split"] for row in rows}, {"development", "test"})


if __name__ == "__main__":
    unittest.main()
