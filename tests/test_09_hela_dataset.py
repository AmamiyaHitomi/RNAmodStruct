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


class HeLaDatasetTests(unittest.TestCase):
    def test_site_table_is_unique_and_combined_ratio_is_count_weighted(self):
        rows = list(read_rows(FINAL / "09_hela_site_level_dataset.csv.gz"))
        self.assertGreater(len(rows), 500)
        self.assertEqual(len({row["site_id"] for row in rows}), len(rows))
        for row in rows:
            expected = int(row["combined_acov"]) / int(row["combined_agcov"])
            self.assertAlmostEqual(float(row["combined_ratio"]), expected, places=7)

    def test_main_dataset_invariants(self):
        rows = list(read_rows(FINAL / "09_hela_main_analysis_dataset.csv.gz"))
        self.assertGreater(len(rows), 500)
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
        rows = read_rows(FINAL / "09_hela_main_analysis_dataset.csv.gz")
        model_ready = 0
        for row in rows:
            if row["model_dataset_included"] == "True":
                model_ready += 1
                self.assertEqual(len(row["sequence_window_201"]), 201)
                self.assertEqual(row["sequence_window_201"][100], "A")
        self.assertGreater(model_ready, 500)

    def test_hela_abundance_is_missing(self):
        rows = read_rows(FINAL / "09_hela_main_analysis_dataset.csv.gz")
        self.assertTrue(all(row["icshape_abundance_rpkm"] == "" for row in rows))

    def test_all_isoform_table_has_one_representative_per_site(self):
        representatives = Counter()
        sites = set()
        for row in read_rows(FINAL / "09_hela_all_isoform_mappings.csv.gz"):
            key = (row["hg38_chr"], row["hg38_pos_1based"], row["hg38_strand"])
            sites.add(key)
            if row["representative_isoform"] == "True":
                representatives[key] += 1
                self.assertEqual(row["isoform_selection_rank"], "1")
        self.assertEqual(set(representatives), sites)
        self.assertTrue(all(count == 1 for count in representatives.values()))

    def test_run_status_passed_without_center_base_failures(self):
        with (ROOT / "results" / "09_hela_build_analysis_dataset_run_status.csv").open(newline="", encoding="utf-8") as handle:
            row = next(csv.DictReader(handle))
        self.assertEqual(row["status"], "PASS")
        self.assertEqual(int(row["center_base_failures"]), 0)


if __name__ == "__main__":
    unittest.main()
