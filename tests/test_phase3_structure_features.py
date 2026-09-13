import gzip
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


class Phase3StructureFeatureTests(unittest.TestCase):
    def test_stage25_gate_and_scale_separation(self):
        status = pd.read_csv(ROOT / "results/status/25_multimodification_structure_features_run_status.csv").iloc[0]
        self.assertEqual(status.status, "PASS")
        self.assertEqual(status.phase3_gate, "GO_TO_STRATIFIED_ASSOCIATION_ANALYSIS")
        self.assertFalse(bool(status.raw_scales_pooled))
        self.assertGreaterEqual(status.distinct_modifications_admitted, 2)

    def test_track_specific_contracts(self):
        attrition = pd.read_csv(ROOT / "results/tables/25_multimodification_structure_attrition.csv")
        self.assertEqual(len(attrition), 5)
        modification_totals = attrition.groupby("modification").main_valid_records.sum()
        self.assertTrue((modification_totals >= 20).all())
        hek_m7g = attrition[attrition.track_id.eq("m7g_hek293t")].iloc[0]
        self.assertEqual(hek_m7g.track_admission, "DESCRIPTIVE_ONLY_LOW_N")
        self.assertTrue((attrition.loc[~attrition.track_id.eq("m7g_hek293t"), "track_admission"] == "INFERENTIAL_READY").all())
        self.assertEqual(set(attrition.measurement_scale), {"quantitative_site_fraction", "semiquantitative_interval_peak", "binary_site"})
        m5c = pd.read_csv(ROOT / "data/final/25_m5c_hela_structure_dataset.csv.gz")
        m7g = pd.read_csv(ROOT / "data/final/25_m7g_hela_structure_dataset.csv.gz")
        self.assertIn("methy_rate_mean", m5c.columns)
        self.assertIn("reactivity_mean_flank10", m5c.columns)
        self.assertIn("fold_enrichment", m7g.columns)
        self.assertIn("reactivity_mean_peak", m7g.columns)
        self.assertNotIn("icshape_center", m7g.columns)

    def test_nm_coordinate_sensitivity(self):
        sensitivity = pd.read_csv(ROOT / "results/tables/25_nm_coordinate_sensitivity.csv")
        self.assertEqual(set(sensitivity.coordinate_offset_nt), {-1, 0, 1})
        self.assertTrue(sensitivity.stability_pass.all())
        self.assertTrue((sensitivity.valid_count_stability_ratio_across_offsets >= 0.90).all())

    def test_all_outputs_are_readable(self):
        for name in [
            "25_m5c_hela_structure_dataset.csv.gz", "25_m7g_hela_structure_dataset.csv.gz",
            "25_m7g_hek293t_structure_dataset.csv.gz", "25_nm_hela_structure_dataset.csv.gz",
            "25_nm_hek293t_structure_dataset.csv.gz",
        ]:
            with gzip.open(ROOT / "data/final" / name, "rt", encoding="utf-8") as handle:
                self.assertTrue(handle.readline().strip())
                self.assertTrue(handle.readline().strip())


if __name__ == "__main__":
    unittest.main()
