import importlib.util
import json
import math
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge


ROOT = Path(__file__).resolve().parents[1]


def load_stage08():
    path = ROOT / "src" / "08_run_prediction_modeling.py"
    spec = importlib.util.spec_from_file_location("stage08_smoke", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def synthetic_frame(rows=30):
    bases = "ACGT"
    records = []
    for index in range(rows):
        prefix = "".join(bases[(index // (4 ** power)) % 4] for power in range(4))
        sequence = (prefix + "ACGT" * 51)[:201]
        record = {
            "site_id": f"site-{index}",
            "analysis_gene": f"gene-{index}",
            "sequence_window_201": sequence,
            "combined_ratio": 0.1 + 0.8 * index / (rows - 1),
            "gc_fraction_21": 0.35 + 0.01 * (index % 10),
            "transcript_position_fraction": (index + 1) / (rows + 1),
            "icshape_abundance_rpkm": 1.0 + index,
            "combined_agcov": 30 + index,
            "coverage_flank10": 0.7 + 0.01 * (index % 10),
            "distance_to_stop_codon_tx": np.nan if index % 7 == 0 else index - 15,
            "distance_to_nearest_splice_edge_tx": np.nan if index % 11 == 0 else index + 3,
            "drach_subtype": "GGACA" if index % 2 else "AAACA",
            "transcript_region": "CDS" if index % 3 else "3UTR",
            "reactivity_mean_flank10": 0.15 + 0.01 * index,
            "predicted_unpaired_probabilities_json": json.dumps([0.2 + 0.01 * (index % 5)] * 201),
        }
        for name in [
            "mfe_per_nt", "ensemble_free_energy_per_nt", "ensemble_diversity",
            "mean_pairing_state_entropy", "predicted_unpaired_mean_up10",
            "predicted_unpaired_mean_down10", "predicted_unpaired_mean_flank10",
            "predicted_unpaired_mean_far",
        ]:
            record[name] = 0.1 + 0.001 * index
        records.append(record)
    return pd.DataFrame(records)


class PipelineSmokeTests(unittest.TestCase):
    def test_group_split_encode_fit_predict_end_to_end(self):
        stage08 = load_stage08()
        frame = synthetic_frame()
        frame["joint_group"] = stage08.make_joint_groups(frame)
        frame["split"] = stage08.assign_split(frame["joint_group"], 0.2, 20260910)
        development = frame[frame["split"].eq("development")]
        test = frame[frame["split"].eq("test")]
        self.assertFalse(set(development["analysis_gene"]) & set(test["analysis_gene"]))
        encoder = stage08.FeatureEncoder().fit(development)
        x_train, names = encoder.transform(development, "M4")
        x_test, test_names = encoder.transform(test, "M4")
        self.assertEqual(names, test_names)
        estimator = Ridge(alpha=300.0).fit(x_train, development["combined_ratio"])
        prediction = np.clip(estimator.predict(x_test), 0.0, 1.0)
        score = stage08.metrics(test["combined_ratio"].to_numpy(), prediction)
        self.assertTrue(all(math.isfinite(value) for value in score.values()))
        self.assertTrue(np.all((prediction >= 0.0) & (prediction <= 1.0)))


if __name__ == "__main__":
    unittest.main()
