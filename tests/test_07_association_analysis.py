import csv
import json
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"


def read_rows(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class AssociationAnalysisTests(unittest.TestCase):
    def test_run_status_and_primary_population(self):
        status = read_rows(ROOT / "results" / "07_hek293t_association_analysis_run_status.csv")[0]
        self.assertEqual(status["status"], "PASS")
        self.assertEqual(int(status["sites"]), 4409)
        self.assertEqual(int(status["genes"]), 1324)
        self.assertEqual(int(status["bootstrap_successful"]), 1000)
        self.assertEqual(int(status["sensitivity_models"]), 7)

    def test_primary_effect_is_finite_and_intervals_exclude_zero(self):
        row = read_rows(TABLES / "07_hek293t_primary_association.csv")[0]
        numeric = [
            "beta_per_sd",
            "ci95_low_per_sd",
            "ci95_high_per_sd",
            "bootstrap_ci95_low_per_sd",
            "bootstrap_ci95_high_per_sd",
        ]
        self.assertTrue(all(math.isfinite(float(row[name])) for name in numeric))
        self.assertGreater(float(row["ci95_low_per_sd"]), 0)
        self.assertGreater(float(row["bootstrap_ci95_low_per_sd"]), 0)
        self.assertLess(float(row["p_value"]), 0.05)

    def test_sensitivity_family_is_complete_and_adjusted(self):
        rows = read_rows(TABLES / "07_hek293t_sensitivity_associations.csv")
        self.assertEqual(len(rows), 7)
        self.assertTrue(all(float(row["beta_per_sd"]) > 0 for row in rows))
        self.assertTrue(all(0 <= float(row["p_adjust_bh_sensitivity_family"]) <= 1 for row in rows))

    def test_profile_has_four_complete_101_nt_grids(self):
        rows = read_rows(TABLES / "07_hek293t_reactivity_profile_by_ratio_quartile.csv")
        self.assertEqual(len(rows), 404)
        for quartile in ["Q1", "Q2", "Q3", "Q4"]:
            positions = sorted(int(row["relative_position"]) for row in rows if row["ratio_quartile"] == quartile)
            self.assertEqual(positions, list(range(-50, 51)))

    def test_posthoc_region_rows_are_explicitly_labelled(self):
        rows = read_rows(TABLES / "07_hek293t_posthoc_region_diagnostic.csv")
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(row["analysis_role"] == "posthoc_fallacy_diagnostic" for row in rows))
        self.assertTrue(all(0 <= float(row["p_adjust_bh_region_diagnostic"]) <= 1 for row in rows))

    def test_figure_qa_has_no_alignment_or_collision_findings(self):
        qa = ROOT / "results" / "figures" / "qa"
        for name in [
            "07_hek293t_association_overview_alignment.json",
            "07_hek293t_association_overview_collision.json",
        ]:
            report = json.loads((qa / name).read_text(encoding="utf-8"))
            self.assertEqual(report["verdict"], "PASS")
            self.assertEqual(int(report["summary"]["fail"]), 0)
            self.assertEqual(int(report["summary"]["warn"]), 0)


if __name__ == "__main__":
    unittest.main()
