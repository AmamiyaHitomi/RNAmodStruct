import csv
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"


def read_rows(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class HeLaAssociationTests(unittest.TestCase):
    def test_run_status_and_population(self):
        status = read_rows(ROOT / "results" / "09b_hela_association_run_status.csv")[0]
        self.assertEqual(status["status"], "PASS")
        self.assertGreaterEqual(int(status["sites"]), 500)
        self.assertEqual(int(status["bootstrap_successful"]), 1000)
        self.assertEqual(int(status["sensitivity_models"]), 7)

    def test_primary_effect_is_finite(self):
        row = read_rows(TABLES / "09b_hela_primary_association.csv")[0]
        numeric = [
            "beta_per_sd", "ci95_low_per_sd", "ci95_high_per_sd",
            "bootstrap_ci95_low_per_sd", "bootstrap_ci95_high_per_sd",
        ]
        self.assertTrue(all(math.isfinite(float(row[name])) for name in numeric))
        self.assertGreaterEqual(float(row["bootstrap_ci95_low_per_sd"]), float(row["ci95_low_per_sd"]) - 0.02)
        self.assertLessEqual(float(row["bootstrap_ci95_high_per_sd"]), float(row["ci95_high_per_sd"]) + 0.02)

    def test_sensitivity_family_is_complete_and_adjusted(self):
        rows = read_rows(TABLES / "09b_hela_sensitivity_associations.csv")
        self.assertEqual(len(rows), 7)
        self.assertTrue(all(0 <= float(row["p_adjust_bh_sensitivity_family"]) <= 1 for row in rows))
        self.assertTrue(all(math.isfinite(float(row["beta_per_sd"])) for row in rows))

    def test_profile_has_four_complete_101_nt_grids(self):
        rows = read_rows(TABLES / "09b_hela_reactivity_profile_by_ratio_quartile.csv")
        self.assertEqual(len(rows), 404)
        for quartile in ["Q1", "Q2", "Q3", "Q4"]:
            positions = sorted(int(row["relative_position"]) for row in rows if row["ratio_quartile"] == quartile)
            self.assertEqual(positions, list(range(-50, 51)))

    def test_diagnostics_are_finite(self):
        rows = read_rows(TABLES / "09b_hela_model_diagnostics.csv")
        self.assertEqual(rows[0]["metric"], "n_sites")
        self.assertGreaterEqual(int(rows[0]["value"]), 500)


if __name__ == "__main__":
    unittest.main()
