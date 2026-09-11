import csv
import gzip
import json
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class PredictedStructureTests(unittest.TestCase):
    def test_run_status(self):
        with (ROOT / "results" / "08a_hek293t_predicted_structure_run_status.csv").open(
            newline="", encoding="utf-8"
        ) as handle:
            row = next(csv.DictReader(handle))
        self.assertEqual(row["status"], "PASS")
        self.assertEqual(int(row["sites"]), 4096)
        self.assertEqual(int(row["unique_sequences"]), 4092)
        self.assertEqual(int(row["probability_vector_length_min"]), 201)
        self.assertEqual(int(row["probability_vector_length_max"]), 201)

    def test_features_are_finite_and_probabilities_are_bounded(self):
        path = ROOT / "data" / "final" / "08a_hek293t_predicted_structure_features.csv.gz"
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 4096)
        self.assertEqual(len({row["site_id"] for row in rows}), 4096)
        for row in rows:
            values = json.loads(row["predicted_unpaired_probabilities_json"])
            self.assertEqual(len(values), 201)
            self.assertTrue(all(math.isfinite(value) and 0 <= value <= 1 for value in values))
            self.assertTrue(math.isfinite(float(row["mfe_per_nt"])))
            self.assertTrue(math.isfinite(float(row["mean_pairing_state_entropy"])))


if __name__ == "__main__":
    unittest.main()
