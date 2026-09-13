import importlib.util
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def load_stage29():
    spec = importlib.util.spec_from_file_location("stage29", ROOT / "src/29_run_phase4_directionality.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Phase4DirectionalityUnitTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.stage29 = load_stage29()

    def test_bh_adjust_is_monotone_in_rank(self):
        adjusted = self.stage29.bh_adjust([0.04, 0.001, 0.03])
        self.assertAlmostEqual(adjusted[1], 0.003)
        self.assertGreaterEqual(adjusted[0], adjusted[2])

    def test_bootstrap_ci_reproducible_and_contains_constant(self):
        values = np.full(30, 0.25)
        first = self.stage29.bootstrap_median_ci(values, 100, 7, 0.05)
        second = self.stage29.bootstrap_median_ci(values, 100, 7, 0.05)
        self.assertEqual(first, second)
        self.assertEqual(first, (0.25, 0.25))

    def test_sign_flip_detects_large_shift(self):
        values = np.linspace(0.5, 1.5, 80)
        p_value = self.stage29.sign_flip_mean_p(values, 1000, 9)
        self.assertLess(p_value, 0.01)


if __name__ == "__main__":
    unittest.main()
