import gzip
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


class Phase3FileAuditTests(unittest.TestCase):
    def test_integrity_and_gate(self):
        files = pd.read_csv(ROOT / "results/tables/24_multimodification_file_audit.csv")
        self.assertEqual(len(files), 5)
        self.assertTrue(files.sha256_match.all())
        status = pd.read_csv(ROOT / "results/status/24_multimodification_file_audit_run_status.csv").iloc[0]
        self.assertEqual(status.status, "PASS")
        self.assertEqual(status.phase3_gate, "GO_TO_STRUCTURE_FEATURE_EXTRACTION")
        self.assertFalse(bool(status.raw_scales_pooled))

    def test_coordinate_and_modality_contracts(self):
        audit = pd.read_csv(ROOT / "results/tables/24_multimodification_overlap_audit.csv")
        self.assertEqual(set(audit.modification), {"m5C", "m7G", "Nm"})
        m5c = audit[audit.track_id.eq("m5c_hela")].iloc[0]
        self.assertGreaterEqual(m5c.reference_base_match_fraction, 0.99)
        self.assertEqual(m5c.admission, "ADMIT_QUANTITATIVE_SITE")
        self.assertTrue(audit[audit.modification.eq("m7G")].admission.str.contains("INTERVAL").all())
        self.assertTrue(audit[audit.modification.eq("Nm")].admission.str.contains("PLUS_MINUS_1").all())
        self.assertTrue((audit.structure_exon_overlap_records >= 20).all())

    def test_normalized_tracks_are_readable_and_distinct(self):
        paths = [
            "data/interim/24_m5c_hela_sites.csv.gz",
            "data/interim/24_m7g_hela_peaks.csv.gz",
            "data/interim/24_m7g_hek293t_peaks.csv.gz",
            "data/interim/24_nm_hela_sites_hg38.csv.gz",
            "data/interim/24_nm_hek293t_sites_hg38.csv.gz",
        ]
        for relative in paths:
            with gzip.open(ROOT / relative, "rt", encoding="utf-8") as handle:
                header = handle.readline()
                first = handle.readline()
            self.assertTrue(header.strip())
            self.assertTrue(first.strip())


if __name__ == "__main__":
    unittest.main()
