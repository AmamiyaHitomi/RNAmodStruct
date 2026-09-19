"""Recompute frozen HEK293T groups and split from processed site records."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "presubmission_stage4" / "split_reproduction.json"
SOURCE = ROOT / "data" / "final" / "06_hek293t_main_analysis_dataset.csv.gz"
ASSIGNMENTS = ROOT / "results" / "tables" / "08_hek293t_group_assignments.csv"
CONFIG = ROOT / "config" / "08_prediction_modeling.yaml"
CODE = ROOT / "src" / "08_run_prediction_modeling.py"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite {OUTPUT}")
    spec = importlib.util.spec_from_file_location("frozen_stage08", CODE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    sites = pd.read_csv(SOURCE, low_memory=False)
    sites = sites.loc[sites.model_dataset_included.astype(str).eq("True")].copy().reset_index(drop=True)
    sites["joint_group"] = module.make_joint_groups(sites)
    sites["split"] = module.assign_split(sites.joint_group, config["split"]["test_fraction"],
                                          config["split"]["random_seed"])
    frozen = pd.read_csv(ASSIGNMENTS, usecols=["site_id", "joint_group", "split"])
    merged = sites[["site_id", "joint_group", "split"]].merge(
        frozen, on="site_id", suffixes=("_replayed", "_frozen"),
        how="outer", indicator=True, validate="one_to_one")
    if not merged._merge.eq("both").all():
        raise AssertionError("Site IDs differ from frozen split")
    for name in ("joint_group", "split"):
        if not merged[f"{name}_replayed"].equals(merged[f"{name}_frozen"]):
            raise AssertionError(f"Frozen {name} not reproduced")
    status = {"status": "PASS", "level": "B_processed_sites_to_group_components_and_split",
              "sites": len(sites), "development": int(sites.split.eq("development").sum()),
              "test": int(sites.split.eq("test").sum()),
              "joint_groups": int(sites.joint_group.nunique()),
              "input_sha256": {str(p.relative_to(ROOT)).replace("\\", "/"): digest(p)
                               for p in (SOURCE, ASSIGNMENTS, CONFIG, CODE)}}
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(status, indent=2), encoding="utf-8")
    print(json.dumps(status))


if __name__ == "__main__":
    main()
