import hashlib
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class Phase4CompletionTests(unittest.TestCase):
    def test_all_stage_statuses_pass(self):
        for stage, name in {
            28: "28_phase4_directionality_run_status.csv",
            29: "29_phase4_directionality_analysis_run_status.csv",
            30: "30_phase4_evidence_and_model_gate_run_status.csv",
            31: "31_phase4_completion_run_status.csv",
        }.items():
            row = pd.read_csv(ROOT / "results" / "status" / name).iloc[0]
            self.assertEqual(int(row.stage), stage)
            self.assertEqual(row.status, "PASS")

    def test_evidence_ceiling_and_model_gate(self):
        status = pd.read_csv(ROOT / "results/status/31_phase4_completion_run_status.csv").iloc[0]
        self.assertEqual(status.highest_evidence_layer, "Directionality")
        self.assertFalse(bool(status.causal_evidence_identified))
        self.assertFalse(bool(status.mechanism_identified))
        self.assertFalse(bool(status.advanced_models_started))

    def test_legacy_candidate_rejected(self):
        audit = pd.read_csv(ROOT / "results/tables/28_phase4_perturbation_source_audit.csv")
        legacy = audit[audit.candidate.str.contains("GSE52662")].iloc[0]
        self.assertEqual(legacy.decision, "REJECT_CAUSAL_CROSS_CELL_LINE")

    def test_release_manifest_hashes(self):
        manifest = pd.read_csv(ROOT / "metadata/manifests/release/phase4_release_manifest.csv")
        self.assertFalse(manifest.path.duplicated().any())
        self.assertGreaterEqual(len(manifest), 20)
        for row in manifest.itertuples(index=False):
            path = ROOT / row.path
            self.assertTrue(path.is_file(), row.path)
            self.assertEqual(path.stat().st_size, row.bytes, row.path)
            self.assertEqual(sha256_file(path), row.sha256, row.path)


if __name__ == "__main__":
    unittest.main()
