import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


class Phase3StratifiedAssociationTests(unittest.TestCase):
    def test_stage26_gate_and_primary_family(self):
        status = pd.read_csv(ROOT / "results/status/26_stratified_multimodification_associations_run_status.csv").iloc[0]
        self.assertEqual(status.status, "PASS")
        self.assertEqual(status.phase3_gate, "GO_TO_PHASE3_SYNTHESIS")
        self.assertEqual(status.primary_tests_fdr_family_size, 3)
        self.assertFalse(bool(status.raw_scales_pooled))
        results = pd.read_csv(ROOT / "results/tables/26_stratified_primary_associations.csv")
        primary = results[results.evidence_tier.eq("primary")]
        self.assertEqual(set(primary.track_id), {"m5c_hela", "m7g_hela", "nm_hela"})
        self.assertTrue(np.isfinite(primary.q_value_primary).all())

    def test_nm_matching_contract(self):
        pairs = pd.read_csv(ROOT / "data/final/26_nm_matched_structure_pairs.csv.gz")
        self.assertTrue((pairs.absolute_distance_nt >= 50).all())
        self.assertFalse(pairs[["track_id", "transcript_id", "control_tx_pos_1based"]].duplicated().any())
        self.assertEqual(set(pairs.matched_base), {"A", "C", "G", "T"})
        qc = pd.read_csv(ROOT / "results/tables/26_nm_matching_qc.csv")
        self.assertTrue((qc.matched_pairs >= 100).all())
        self.assertTrue((qc.matching_fraction > 0.80).all())

    def test_low_n_track_is_descriptive_only(self):
        low = pd.read_csv(ROOT / "results/tables/26_descriptive_low_n_tracks.csv").iloc[0]
        self.assertEqual(low.track_id, "m7g_hek293t")
        self.assertEqual(low.analysis_role, "DESCRIPTIVE_ONLY_LOW_N")
        self.assertFalse(bool(low.inferential_p_value_reported))

    def test_m5c_replicates_remain_separate(self):
        table = pd.read_csv(ROOT / "results/tables/26_m5c_replicate_sensitivity.csv")
        self.assertEqual(set(table.outcome), {"methy_rate_B", "methy_rate_C", "methy_rate_E"})
        self.assertEqual(len(table), 3)


if __name__ == "__main__":
    unittest.main()
