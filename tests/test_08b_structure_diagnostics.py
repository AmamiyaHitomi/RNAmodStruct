import csv
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"


def read_rows(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class StructureDiagnosticTests(unittest.TestCase):
    def test_run_status_is_development_only(self):
        row = read_rows(ROOT / "results" / "08b_hek293t_structure_diagnostics_run_status.csv")[0]
        self.assertEqual(row["status"], "PASS")
        self.assertEqual(row["analysis_scope"], "development_only")
        self.assertEqual(int(row["sites"]), 3271)
        self.assertEqual(int(row["blocks_compared"]), 5)
        self.assertEqual(row["selected_block"], "R_none")

    def test_all_blocks_have_selected_finite_cv_results(self):
        rows = read_rows(TABLES / "08b_hek293t_expanded_structure_summary.csv")
        self.assertEqual({row["block"] for row in rows}, {
            "R_none", "R_core", "R_directional", "R_profile_no_center", "R_profile_all"
        })
        self.assertTrue(all(row["analysis_scope"] == "development_GroupKFold_only" for row in rows))
        self.assertTrue(all(math.isfinite(float(row["mae"])) for row in rows))

    def test_replicate_agreement_is_valid(self):
        rows = read_rows(TABLES / "08b_hek293t_glori_replicate_agreement.csv")
        overall = next(row for row in rows if row["subset"] == "all_development")
        self.assertEqual(int(overall["sites"]), 3271)
        self.assertGreaterEqual(float(overall["replicate_mae"]), 0)
        self.assertLessEqual(abs(float(overall["replicate_spearman"])), 1)


if __name__ == "__main__":
    unittest.main()
