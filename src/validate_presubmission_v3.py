"""Validate the main_v3 evidence links and A/B presubmission replay."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pymupdf


ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "docs" / "manuscript_v3"
OUT = ROOT / "results" / "presubmission_stage4" / "evidence_validation_v3.json"
checks = []


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def check(name: str, ok: bool, detail: str = "") -> None:
    checks.append({"name": name, "pass": bool(ok), "detail": detail})
    if not ok:
        raise AssertionError(f"{name}: {detail}")


def main() -> None:
    if OUT.exists():
        raise FileExistsError(f"Refusing to overwrite {OUT}")
    main = (PAPER / "main_v3.tex").read_text(encoding="utf-8")
    supplement = (PAPER / "supplement_v3.tex").read_text(encoding="utf-8")
    bibliography = (PAPER / "references.bib").read_text(encoding="utf-8")
    keys = set(re.findall(r"@\w+\{([^,]+)", bibliography))
    cited = set(k.strip() for group in re.findall(r"\\cite\w*\{([^}]+)\}", main) for k in group.split(","))
    check("Main bibliography keys resolve", cited <= keys, str(sorted(cited - keys)))
    check("SMART-m6A continuous regression acknowledged", "GLORI-labeled regression of continuous modification intensity" in main)
    check("Fig 4 conclusion narrowed", "did not improve the evaluated frozen transfer models" in main)
    check("Development-only exclusion disclosed", "does not remove all previously viewed HEK293T sites" in main)
    check("Post hoc matching disclosed", "Matched post hoc comparisons" in main and "post hoc sensitivities" in supplement)
    check("Normalized replicate labels disclosed", "replicate-specific" in supplement and "no simple replicate mean" in supplement)
    check("Stage-19 distinct estimands disclosed", "not a common-raw-unit slope difference" in supplement)
    check("No obsolete broad Fig 4 claim", "consistently reject an experimental-structure gain" not in main)

    manifest = json.loads((PAPER / "source_data_manifest.json").read_text(encoding="utf-8"))
    for filename, expected in manifest.items():
        copied, original = PAPER / "source_data" / filename, ROOT / "results" / "tables" / filename
        check("Frozen source hash: " + filename,
              copied.exists() and original.exists() and sha256(copied) == expected == sha256(original))
    check("Four rendered figure PDFs", all((PAPER / "figures" / f"figure{i}_{name}.pdf").exists()
                                          for i, name in ((1, "attrition"), (2, "association"),
                                                          (3, "internal"), (4, "transfer"))))
    check("Five rebuilt LaTeX tables", len(list((PAPER / "tables").glob("*.tex"))) == 5)

    reproduction = pd.read_csv(ROOT / "results" / "presubmission_stage1" / "published_metric_reproduction.csv")
    check("60 frozen prediction metrics reproduced", len(reproduction) == 60 and
          reproduction.absolute_difference.max() < 1e-8)
    original = pd.read_csv(ROOT / "results" / "tables" / "09b_hela_primary_association.csv").iloc[0]
    sensitivity = pd.read_csv(ROOT / "results" / "presubmission_stage3" / "hela_outcome_sensitivity.csv")
    check("HeLa primary coefficient reproduced", np.isclose(sensitivity.iloc[0].beta_raw, original.beta_raw, atol=1e-12))
    check("Two normalized-label sensitivities retained", set(sensitivity.outcome) ==
          {"combined_ratio", "norme_ratio_rep1", "norme_ratio_rep2"})
    cohorts = json.loads((ROOT / "metadata" / "manifests" / "presubmission_phase0" / "freeze_manifest.json").read_text())
    check("Three frozen HeLa cohorts", [cohorts["cohorts"][name]["site_count"] for name in
                                         ("all_qualifying", "development_overlap_excluded", "hek_main_4409_overlap_excluded")]
          == [24960, 22256, 21351])
    split = json.loads((ROOT / "results" / "presubmission_stage4" / "split_reproduction.json").read_text())
    replay = json.loads((ROOT / "results" / "presubmission_stage4" / "replay" / "replay_validation.json").read_text())
    check("B split replay", split["status"] == "PASS" and split["development"] == 3271 and split["test"] == 825)
    check("B model replay", replay["status"] == "PASS" and replay["saved_model_count"] == 16 and
          len(replay["tables"]) == 6 and max(row["max_abs_numeric_difference"] for row in replay["tables"]) < 1e-8)

    pdf_info = {}
    for stem, minimum in (("main_v3", 8), ("supplement_v3", 4)):
        doc = pymupdf.open(PAPER / (stem + ".pdf"))
        page_text = "\n".join(page.get_text() for page in doc)
        outside = []
        for index, page in enumerate(doc):
            for block in page.get_text("blocks"):
                if len(block) < 5 or not str(block[4]).strip():
                    continue
                x0, y0, x1, y1 = block[:4]
                if x0 < -1 or y0 < -1 or x1 > page.rect.width + 1 or y1 > page.rect.height + 1:
                    outside.append(index + 1)
        check(stem + " page count and text bounds", len(doc) >= minimum and not outside,
              f"pages={len(doc)}, outside={outside}")
        check(stem + " no unresolved references", "??" not in page_text)
        pdf_info[stem] = {"pages": len(doc), "sha256": sha256(PAPER / (stem + ".pdf"))}
    for stem in ("main_v3", "supplement_v3"):
        log = (PAPER / (stem + ".log")).read_text(errors="replace")
        check(stem + " no LaTeX error or overfull box", not re.search(r"(^! |Overfull \\hbox|Undefined control sequence|undefined references)", log, re.I | re.M))

    result = {"status": "PASS", "level_A": "frozen_tables_to_figures_tables_and_PDF",
              "level_B": "processed_sites_to_split_then_processed_input_to_models_and_metrics",
              "check_count": len(checks), "checks": checks, "pdf": pdf_info}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"status": "PASS", "checks": len(checks), "pdf": pdf_info}))


if __name__ == "__main__":
    main()
