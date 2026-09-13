import unittest
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]


class Phase3ScoutTests(unittest.TestCase):
    def test_config_has_one_card_per_required_modification(self):
        config = yaml.safe_load((ROOT / "config/23_multimodification_source_scout.yaml").read_text(encoding="utf-8"))
        cards = config["sources"]
        self.assertEqual(config["stage"], 23)
        self.assertEqual(config["phase"], 3)
        self.assertEqual({row["modification"] for row in cards}, set(config["required_modifications"]))
        self.assertEqual(len(cards), len(config["required_modifications"]))

    def test_scout_gate_and_download_queue(self):
        status = pd.read_csv(ROOT / "results/status/23_multimodification_source_scout_run_status.csv").iloc[0]
        self.assertEqual(status.status, "PASS")
        self.assertEqual(status.phase3_gate, "GO_TO_DOWNLOAD_AUDIT")
        self.assertGreaterEqual(status.download_ready_modifications, 2)
        queue = pd.read_csv(ROOT / "metadata/manifests/download/phase3_download_queue.csv")
        self.assertEqual(len(queue), 5)
        self.assertEqual(set(queue.local_status), {"DOWNLOADED_VERIFIED"})
        self.assertTrue(queue.url.str.startswith("https://ftp.ncbi.nlm.nih.gov/").all())
        self.assertFalse(queue.local_file.duplicated().any())

    def test_measurement_tracks_remain_separate(self):
        cards = pd.read_csv(ROOT / "results/tables/23_multimodification_source_cards.csv")
        decisions = dict(zip(cards.modification, cards.decision))
        self.assertEqual(decisions["m1A"], "HOLD_TECHNICAL_VALIDITY")
        self.assertEqual(decisions["A-to-I"], "HOLD_NO_MATCHED_CELL")
        self.assertIn("SEMIQUANTITATIVE", decisions["m7G"])
        self.assertIn("BINARY", decisions["Nm"])


if __name__ == "__main__":
    unittest.main()
