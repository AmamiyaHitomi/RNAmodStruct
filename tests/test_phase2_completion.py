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


class Phase2CompletionTests(unittest.TestCase):
    def test_required_stages_pass(self):
        patterns = {
            18: "18_public_data_audit_run_status.csv",
            19: "19_invivo_invitro_run_status.csv",
            20: "20_sacseq_cross_technology_run_status.csv",
            21: "21_conditional_source_admission_run_status.csv",
            22: "22_phase2_completion_run_status.csv",
        }
        for stage, name in patterns.items():
            row = pd.read_csv(ROOT / "results" / "status" / name).iloc[0]
            self.assertEqual(int(row.stage), stage)
            self.assertEqual(row.status, "PASS")

    def test_conditional_sources_remain_held(self):
        decisions = pd.read_csv(ROOT / "results/tables/21_conditional_source_admission.csv")
        self.assertEqual(set(decisions.record_id), {"human_rna_map", "pars_gse50676"})
        self.assertEqual(set(decisions.decision), {"HOLD_NOT_ADMITTED"})
        self.assertFalse(decisions.matched_quantitative_m6a_available.any())

    def test_phase2_exit_gates_pass(self):
        gates = pd.read_csv(ROOT / "results/tables/22_phase2_exit_gates.csv")
        self.assertEqual(len(gates), 4)
        self.assertTrue(gates.required.all())
        self.assertTrue(gates.observed.all())
        self.assertEqual(set(gates.status), {"PASS"})

    def test_phase2_manifest_hashes(self):
        manifest = pd.read_csv(ROOT / "metadata/manifests/release/phase2_release_manifest.csv")
        self.assertFalse(manifest.path.duplicated().any())
        self.assertGreaterEqual(len(manifest), 20)
        for row in manifest.itertuples(index=False):
            path = ROOT / row.path
            self.assertTrue(path.is_file(), row.path)
            self.assertEqual(path.stat().st_size, row.bytes, row.path)
            self.assertEqual(sha256_file(path), row.sha256, row.path)


if __name__ == "__main__":
    unittest.main()
