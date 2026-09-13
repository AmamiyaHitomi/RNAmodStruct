import importlib.util
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
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


class Phase1ExtensionTests(unittest.TestCase):
    def test_all_stage_contracts_and_statuses_exist(self):
        for stage in range(12, 18):
            self.assertEqual(len(list((ROOT / "config").glob(f"{stage}_*.yaml"))), 1)
            self.assertEqual(len(list((ROOT / "src").glob(f"{stage}_*.py"))), 1)
            config = yaml.safe_load(next((ROOT / "config").glob(f"{stage}_*.yaml")).read_text(encoding="utf-8"))
            self.assertEqual(int(config["stage"]), stage)
            status = pd.read_csv(next((ROOT / "results" / "status").glob(f"{stage}_*run_status.csv")))
            self.assertEqual(status.loc[0, "status"], "PASS")
            self.assertEqual(len(list((ROOT / "results" / "reports").glob(f"{stage}_*report.md"))), 1)
            self.assertGreaterEqual(len(list((ROOT / "results/tables").glob(f"{stage}_*"))), 1)

    def test_meta_analysis_preserves_dataset_effects_and_two_models(self):
        rows = pd.read_csv(ROOT / "results/tables/12_meta_analysis.csv")
        self.assertEqual(set(rows.loc[rows.row_type.eq("dataset"), "dataset"]), {"HEK293T", "HeLa"})
        self.assertEqual(set(rows.loc[rows.row_type.eq("pooled"), "model"]), {"fixed_effect", "random_effect_DL"})
        pooled = rows[rows.row_type.eq("pooled")]
        self.assertTrue((pooled.ci95_low <= pooled.effect_per_sd).all())
        self.assertTrue((pooled.effect_per_sd <= pooled.ci95_high).all())

    def test_conditional_cross_fitting_has_no_group_fold_leakage(self):
        residuals = pd.read_csv(ROOT / "results/tables/13_cross_fitted_residuals.csv.gz")
        self.assertTrue(np.isfinite(residuals[["residual_m", "residual_r"]]).all().all())
        folds_per_group = residuals.groupby(["dataset", "joint_group"])["fold"].nunique()
        self.assertEqual(int(folds_per_group.max()), 1)
        results = pd.read_csv(ROOT / "results/tables/13_conditional_dependence.csv")
        self.assertTrue(results["group_sign_flip_p_value"].between(0, 1).all())

    def test_all_missing_training_covariate_is_finite(self):
        common = load_script("phase1_common.py")
        rows = 8
        frame = pd.DataFrame({
            "gc_fraction_21": np.linspace(.3, .6, rows), "transcript_position_fraction": np.linspace(.1, .8, rows),
            "icshape_abundance_rpkm": [np.nan] * rows, "combined_agcov": np.arange(rows) + 20,
            "coverage_flank10": 1.0, "distance_to_stop_codon_tx": [np.nan] * rows,
            "distance_to_nearest_splice_edge_tx": np.arange(rows), "drach_subtype": "GGACT",
            "transcript_region": "CDS", "sequence_window_201": ["A" * 201] * rows,
            "reactivity_mean_flank10": np.linspace(.1, .2, rows),
        })
        encoded, _ = common.TabularEncoder().fit(frame).transform(frame, "M1")
        self.assertTrue(np.isfinite(encoded).all())

    def test_entropy_family_and_bh_contract(self):
        rows = pd.read_csv(ROOT / "results/tables/14_structure_entropy_associations.csv")
        self.assertEqual(len(rows), 2 * 3 * 8)
        self.assertEqual(rows.groupby(["dataset", "model"])["predictor"].nunique().min(), 8)
        self.assertTrue(rows.p_value_bh.between(0, 1).all())
        self.assertEqual(int(rows.is_primary_entropy_hypothesis.sum()), 6)

    def test_isoform_estimands_are_site_level(self):
        rows = pd.read_csv(ROOT / "results/tables/15_isoform_sensitivity_associations.csv")
        self.assertEqual(len(rows), 8)
        self.assertTrue((rows.groupby("dataset")["estimand"].nunique() == 4).all())
        sites = pd.read_csv(ROOT / "results/tables/15_isoform_site_estimates.csv.gz")
        self.assertFalse(sites.duplicated(["dataset", "estimand", "site_id"]).any())

    def test_401nt_outputs_and_cache_are_complete(self):
        cache = pd.read_csv(ROOT / "data/interim/16_401nt_vienna_cache.csv")
        self.assertEqual(len(cache), cache.sequence_sha256_401.nunique())
        self.assertFalse(cache.isna().any().any())
        for dataset in ("hek293t", "hela"):
            frame = pd.read_csv(ROOT / f"data/final/16_{dataset}_401nt_dataset.csv.gz")
            self.assertTrue(frame.sequence_window_401.str.len().eq(401).all())
            self.assertTrue(frame.coverage_up200.ge(.70).all() and frame.coverage_down200.ge(.70).all())
        metrics = pd.read_csv(ROOT / "results/tables/16_window_401_external_metrics.csv")
        self.assertEqual(set(metrics.model), {"M0", "M1", "M2", "M3", "M4"})

    def test_401_checkpoint_resumes_without_refolding(self):
        stage16 = load_script("16_run_window_401_sensitivity.py")
        sequences = [("ACGT" * 101)[:400] + base for base in "ACGT"]
        frame = pd.DataFrame({"sequence_window_401": sequences})
        frame["sequence_sha256_401"] = frame["sequence_window_401"].map(
            lambda sequence: stage16.hashlib.sha256(sequence.encode("ascii")).hexdigest()
        )
        config = {"predicted_structure": {"temperature_celsius": 37.0, "dangles": 2}}
        def fake_fold(sequence, temperature, dangles):
            digest = stage16.hashlib.sha256(sequence.encode("ascii")).hexdigest()
            return {"sequence_sha256_401": digest, "mfe_per_nt_401": -.1,
                    "ensemble_free_energy_per_nt_401": -.2, "ensemble_diversity_401": 1.0,
                    "mean_pairing_state_entropy_401": .3, "predicted_unpaired_mean_flank200": .5}
        with tempfile.TemporaryDirectory() as directory:
            cache = Path(directory) / "cache.csv"
            with patch.object(stage16, "CACHE", cache), patch.object(stage16, "fold_401", side_effect=fake_fold) as mocked:
                first = stage16.attach_vienna([frame], config)[0]
                self.assertEqual(mocked.call_count, 4)
                second = stage16.attach_vienna([frame], config)[0]
                self.assertEqual(mocked.call_count, 4)
                self.assertEqual(len(first), len(second))

    def test_nonlinear_matrix_and_explanations(self):
        external = pd.read_csv(ROOT / "results/tables/17_nonlinear_external_metrics.csv")
        self.assertEqual(len(external), 10)
        self.assertEqual(set(external.model), {"M0", "M1", "M2", "M3", "M4"})
        self.assertEqual(set(external.learner), {"Ridge", "HistGradientBoostingRegressor"})
        self.assertTrue(np.isfinite(external[["mae", "rmse", "spearman", "r2"]]).all().all())
        importance = pd.read_csv(ROOT / "results/tables/17_grouped_permutation_importance.csv")
        self.assertEqual(set(importance.feature_block), {"common_covariates", "sequence_3mer", "predicted_structure", "experimental_structure"})
        ablation = pd.read_csv(ROOT / "results/tables/17_retraining_ablation.csv")
        self.assertEqual(len(ablation), 8)

    def test_phase_completion_reproducibility_and_manifest(self):
        reproducibility = pd.read_csv(ROOT / "results/status/phase1_reproducibility_status.csv").iloc[0]
        self.assertEqual(reproducibility["status"], "PASS")
        self.assertEqual(reproducibility["verdict"], "REPRODUCIBLE")
        checks = pd.read_csv(ROOT / "results/tables/phase1_reproducibility_check.csv")
        self.assertEqual(len(checks), 16)
        self.assertEqual(set(checks.verdict), {"MATCH"})
        validation = pd.read_csv(ROOT / "results/tables/phase1_statistical_validation.csv")
        self.assertEqual(len(validation), 11)
        completion = pd.read_csv(ROOT / "results/status/phase1_completion_status.csv").iloc[0]
        self.assertEqual(completion["status"], "PASS")
        self.assertEqual(int(completion["stages_passed"]), 6)
        manifest = pd.read_csv(ROOT / "metadata/manifests/release/phase1_release_manifest.csv")
        self.assertFalse(manifest.path.duplicated().any())
        self.assertTrue(set(f"config/{stage}_" for stage in range(12, 18)).issubset(
            {value[:10] for value in manifest.path if value.startswith("config/")}
        ))
        self.assertNotIn("metadata/manifests/release/release_manifest.csv", set(manifest.path))
        import hashlib
        for row in manifest.itertuples():
            path = ROOT / row.path
            self.assertTrue(path.is_file(), row.path)
            self.assertEqual(path.stat().st_size, row.bytes, row.path)
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), row.sha256, row.path)
        report = (ROOT / "results/reports/phase1_summary_report.md").read_text(encoding="utf-8")
        self.assertIn("Verification Status: VERIFIED", report)


if __name__ == "__main__":
    unittest.main()
