import importlib.util
import sys
import unittest
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))


def load_script(name):
    path = ROOT / "src" / name
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Phase2ExtensionTests(unittest.TestCase):
    def test_stage_configs_are_isolated_from_frozen_mainline(self):
        for stage in (18, 19):
            config = yaml.safe_load(next((ROOT / "config").glob(f"{stage}_*.yaml")).read_text(encoding="utf-8"))
            self.assertEqual(config["stage"], stage)
            self.assertEqual(config["phase"], 2)

    def test_stage18_audit_contract(self):
        rows = pd.read_csv(ROOT / "results/tables/18_public_data_audit.csv")
        local = rows[rows.local_file.notna() & rows.local_file.ne("")]
        self.assertEqual(len(local), 3)
        self.assertEqual(set(local.integrity), {"GZIP_FULL_STREAM_OK"})
        self.assertEqual(len(local.sha256.str.fullmatch(r"[0-9a-f]{64}")), 3)
        bed = local[local.local_file.str.endswith(".bed.gz")]
        self.assertTrue(bed.coordinate_checks.str.contains("width1=True;strand=True;unique=True", regex=False).all())

    def test_flank_metrics_exclude_center(self):
        stage19 = load_script("19_run_invivo_invitro_pairing.py")
        values = stage19.array("f", range(1, 22))
        metrics = stage19.flank_metrics(values, 11, 10)
        self.assertAlmostEqual(metrics["coverage_flank10"], 1.0)
        self.assertAlmostEqual(metrics["reactivity_mean_flank10"], sum(list(range(1, 11)) + list(range(12, 22))) / 20)

    def test_stage19_population_and_effect_contract(self):
        paired = pd.read_csv(ROOT / "data/final/19_hek293t_invivo_invitro_paired.csv.gz")
        self.assertFalse(paired.site_id.duplicated().any())
        self.assertTrue(paired.coverage_up10_invivo.ge(.70).all())
        self.assertTrue(paired.coverage_down10_invivo.ge(.70).all())
        self.assertTrue(paired.coverage_up10_invitro.ge(.70).all())
        self.assertTrue(paired.coverage_down10_invitro.ge(.70).all())
        effects = pd.read_csv(ROOT / "results/tables/19_paired_condition_associations.csv")
        self.assertEqual(set(effects.row_type), {"condition_effect", "paired_effect_contrast"})
        self.assertEqual(set(effects.loc[effects.row_type.eq("condition_effect"), "analysis"]), {"in_vivo", "in_vitro", "delta_invivo_minus_invitro"})


if __name__ == "__main__":
    unittest.main()
