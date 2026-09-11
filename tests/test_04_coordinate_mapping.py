import csv
import math
import sys
import unittest
from array import array
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pipeline_common import ChainLiftOver, Exon, coverage, genome_to_tx, tx_to_genome


class CoordinateMappingTests(unittest.TestCase):
    def test_positive_strand_across_exons(self):
        exons = [
            Exon("chr1", 100, 109, "+", "tx", "g", "G", 1, 1, 10),
            Exon("chr1", 200, 209, "+", "tx", "g", "G", 2, 11, 20),
        ]
        self.assertEqual(tx_to_genome(exons, 12), ("chr1", 201, "+", 2))
        self.assertEqual(genome_to_tx(exons[1], 201), 12)

    def test_negative_strand_across_exons(self):
        exons = [
            Exon("chr1", 200, 209, "-", "tx", "g", "G", 1, 1, 10),
            Exon("chr1", 100, 109, "-", "tx", "g", "G", 2, 11, 20),
        ]
        self.assertEqual(tx_to_genome(exons, 12), ("chr1", 108, "-", 2))
        self.assertEqual(genome_to_tx(exons[1], 108), 12)

    def test_window_coverage_excludes_center(self):
        values = array("f", [1.0] * 21)
        values[10] = math.nan
        left, right, combined = coverage(values, 11, 10)
        self.assertEqual((left, right, combined), (1.0, 1.0, 1.0))

    def test_recorded_chain_mapping_roundtrip(self):
        checks = ROOT / "metadata" / "coordinate" / "04_hek293t_coordinate_checks.csv"
        with checks.open(newline="", encoding="utf-8") as handle:
            row = next(r for r in csv.DictReader(handle) if r["mapping_status"] == "unique")
        lifter = ChainLiftOver(ROOT / "data" / "reference" / "ucsc" / "hg19ToHg38.over.chain.gz")
        forward = lifter.forward(row["hg19_chr"], int(row["hg19_pos_1based"]), row["hg19_strand"])
        self.assertTrue(any(c == row["hg38_chr"] and p == int(row["hg38_pos_1based"]) and s == row["hg38_strand"] for c, p, s, _ in forward))
        reverse = lifter.reverse(row["hg38_chr"], int(row["hg38_pos_1based"]), row["hg38_strand"])
        self.assertTrue(any(c == row["hg19_chr"] and p == int(row["hg19_pos_1based"]) and s == row["hg19_strand"] for c, p, s, _ in reverse))


if __name__ == "__main__":
    unittest.main()
