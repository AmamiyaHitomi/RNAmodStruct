import csv
import gzip
import json
import unittest
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FINAL = ROOT / "data" / "final"


def read_rows(path):
    with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle)


class AnalysisDatasetTests(unittest.TestCase):
    def test_site_table_is_unique_and_combined_ratio_is_count_weighted(self):
        rows = list(read_rows(FINAL / "06_hek293t_site_level_dataset.csv.gz"))
        self.assertEqual(len(rows), 20559)
        self.assertEqual(len({row["site_id"] for row in rows}), len(rows))
        for row in rows:
            expected = int(row["combined_acov"]) / int(row["combined_agcov"])
            self.assertAlmostEqual(float(row["combined_ratio"]), expected, places=7)

    def test_main_dataset_invariants(self):
        rows = list(read_rows(FINAL / "06_hek293t_main_analysis_dataset.csv.gz"))
        self.assertEqual(len(rows), 4409)
        for row in rows:
            self.assertEqual(row["detected_both_replicates"], "True")
            self.assertEqual(row["main_window_valid"], "True")
            self.assertEqual(row["main_analysis_included"], "True")
            self.assertGreaterEqual(int(row["valid_count_up10"]), 7)
            self.assertGreaterEqual(int(row["valid_count_down10"]), 7)
            window = json.loads(row["reactivity_window_m50_p50_json"])
            mask = json.loads(row["reactivity_mask_m50_p50_json"])
            self.assertEqual(len(window), 101)
            self.assertEqual(len(mask), 101)
            self.assertEqual(sum(mask), int(row["valid_reactivity_count_m50_p50"]))

    def test_full_sequence_windows_have_center_a(self):
        rows = read_rows(FINAL / "06_hek293t_main_analysis_dataset.csv.gz")
        model_ready = 0
        for row in rows:
            if row["model_dataset_included"] == "True":
                model_ready += 1
                self.assertEqual(len(row["sequence_window_201"]), 201)
                self.assertEqual(row["sequence_window_201"][100], "A")
        self.assertEqual(model_ready, 4096)

    def test_all_isoform_table_has_one_representative_per_site(self):
        representatives = Counter()
        sites = set()
        total = 0
        for row in read_rows(FINAL / "06_hek293t_all_isoform_mappings.csv.gz"):
            total += 1
            key = (row["hg38_chr"], row["hg38_pos_1based"], row["hg38_strand"])
            sites.add(key)
            if row["representative_isoform"] == "True":
                representatives[key] += 1
                self.assertEqual(row["isoform_selection_rank"], "1")
        self.assertEqual(total, 33962)
        self.assertEqual(len(sites), 20559)
        self.assertEqual(set(representatives), sites)
        self.assertTrue(all(count == 1 for count in representatives.values()))


if __name__ == "__main__":
    unittest.main()
