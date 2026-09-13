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


class Phase3CompletionTests(unittest.TestCase):
    def test_required_stages_and_completion_pass(self):
        names = {
            23: "23_multimodification_source_scout_run_status.csv",
            24: "24_multimodification_file_audit_run_status.csv",
            25: "25_multimodification_structure_features_run_status.csv",
            26: "26_stratified_multimodification_associations_run_status.csv",
            27: "27_phase3_completion_run_status.csv",
        }
        for stage, name in names.items():
            row = pd.read_csv(ROOT / "results" / "status" / name).iloc[0]
            self.assertEqual(int(row.stage), stage)
            self.assertEqual(row.status, "PASS")

    def test_exit_gates_and_next_phase(self):
        gates = pd.read_csv(ROOT / "results/tables/27_phase3_exit_gates.csv")
        self.assertEqual(len(gates), 6)
        self.assertTrue(gates.required.all())
        self.assertTrue(gates.observed.all())
        status = pd.read_csv(ROOT / "results/status/27_phase3_completion_run_status.csv").iloc[0]
        self.assertEqual(status.next_phase_gate, "PHASE3_COMPLETE_GO_TO_PHASE4_PLANNING")
        self.assertEqual(status.primary_findings_fdr_significant, 0)

    def test_claim_boundaries_prevent_overstatement(self):
        claims = pd.read_csv(ROOT / "results/tables/27_phase3_claim_boundaries.csv")
        self.assertEqual(len(claims), 5)
        self.assertFalse(claims.causal_claim_permitted.any())
        self.assertFalse(claims.exact_zero_claim_permitted.any())
        self.assertTrue(claims.permitted_statement.str.len().gt(20).all())
        self.assertTrue(claims.forbidden_statement.str.len().gt(20).all())

    def test_release_manifest_hashes(self):
        manifest = pd.read_csv(ROOT / "metadata/manifests/release/phase3_release_manifest.csv")
        self.assertFalse(manifest.path.duplicated().any())
        self.assertGreaterEqual(len(manifest), 35)
        for row in manifest.itertuples(index=False):
            path = ROOT / row.path
            self.assertTrue(path.is_file(), row.path)
            self.assertEqual(path.stat().st_size, row.bytes, row.path)
            self.assertEqual(sha256_file(path), row.sha256, row.path)


if __name__ == "__main__":
    unittest.main()
