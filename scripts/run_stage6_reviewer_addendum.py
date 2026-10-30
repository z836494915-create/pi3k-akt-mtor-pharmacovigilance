from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
from datetime import datetime
from pathlib import Path

from statsmodels.stats.multitest import multipletests


from package_paths import OUTPUT

ROOT = OUTPUT
STAGE5 = ROOT / "03_确认后补充分析" / "阶段五_主要分析"
OUT = ROOT / "03_确认后补充分析" / "阶段六_敏感性分析"
STAGE5_SCRIPT = Path(__file__).with_name("run_stage5_primary_analysis.py")

RENAL_PTS = {
    "PT_acute_kidney_injury": "Acute kidney injury",
    "PT_renal_failure": "Renal failure",
    "PT_renal_impairment": "Renal impairment",
}
DEHYDRATION_PT = "PT_dehydration"


def load_stage5_module():
    spec = importlib.util.spec_from_file_location("stage5_primary_addendum", STAGE5_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load Stage 5 definitions from {STAGE5_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.OUTCOMES[DEHYDRATION_PT] = ["Dehydration"]
    for key, pt in RENAL_PTS.items():
        module.OUTCOMES[key] = [pt]
    return module


def write_csv(path: Path, rows: list[dict]):
    if not rows:
        raise RuntimeError(f"Refusing to write an empty results file: {path}")
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def apply_by(rows: list[dict], q_name: str):
    # Stage 5 now includes the renal PTs. Replace those hypotheses by the
    # current renal results instead of counting them a second time.
    with (STAGE5 / "阶段五_次要结局ROR_Fisher_BY.csv").open(encoding="utf-8-sig", newline="") as f:
        prior = list(csv.DictReader(f))
    key = lambda r: (r["database"], r["drug"], r["outcome"])
    if len({key(r) for r in prior}) != len(prior):
        raise RuntimeError("Duplicate hypothesis keys in Stage 5 secondary results")
    if len({key(r) for r in rows}) != len(rows):
        raise RuntimeError("Duplicate hypothesis keys in current renal results")
    if any(r["outcome"] not in RENAL_PTS for r in rows):
        raise RuntimeError("Unexpected outcome in renal PT results")
    prior_eligible = [
        r
        for r in prior
        if r["outcome"] not in RENAL_PTS and r["estimable"].strip().lower() == "true"
    ]
    current_eligible = [r for r in rows if str(r["estimable"]).strip().lower() == "true"]
    if len(prior_eligible) != 98 or len(current_eligible) != 42 or len(rows) != 45:
        raise RuntimeError(
            f"Unexpected secondary family: non-renal={len(prior_eligible)}, "
            f"renal={len(current_eligible)}/{len(rows)}; expected 98+42=140"
        )
    if len({key(r) for r in prior_eligible + current_eligible}) != 140:
        raise RuntimeError("Combined secondary family must contain 140 unique estimable hypotheses")
    all_p = [float(r["fisher_two_sided_p"]) for r in prior_eligible]
    all_p.extend(float(r["_p"]) for r in current_eligible)
    if all_p:
        q_values = multipletests(all_p, method="fdr_by")[1]
        q_current = q_values[len(prior_eligible) :]
        for row, q in zip(current_eligible, q_current):
            row[q_name] = f"{float(q):.6g}"
    return len(prior_eligible), len(prior_eligible) + len(current_eligible)


def public_rows(rows: list[dict]) -> list[dict]:
    return [{k: v for k, v in row.items() if not k.startswith("_")} for row in rows]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    s5 = load_stage5_module()
    loaded = []
    for loader in (s5.load_faers, s5.load_canada, s5.load_jader):
        print(f"Loading {loader.__name__}...", flush=True)
        loaded.append(loader())

    stage5_primary = {}
    with (STAGE5 / "阶段五_主要结局ROR_Fisher_Holm.csv").open(encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            stage5_primary[(row["database"], row["drug"])] = row

    renal_rows = []
    dka_dehydration_rows = []
    renal_overlap_rows = []
    reconciliation = []
    source_records = []

    for n_reports, targets, events, meta in loaded:
        db = meta["database"]
        primary = events["primary_union"]
        dka = events["PT_diabetic_ketoacidosis"]
        strict_dehydration = events[DEHYDRATION_PT]
        dehydration_domain = events["dehydration_electrolyte_domain"]

        for source in meta.get("source_zips", meta.get("source_files", [])):
            base = (
                s5.FAERS_ZIPS
                if db == "FAERS"
                else (s5.CAN if db == "Canada Vigilance" else s5.JADER)
            )
            source_path = base / source["name"]
            source_records.append({"database": db, **source, "path": str(source_path)})

        for drug, exposed in targets.items():
            locked = stage5_primary[(db, drug)]
            loaded_target_n = len(exposed)
            loaded_primary_n = len(exposed & primary)
            ok = (
                loaded_target_n == int(locked["target_reports_n"])
                and loaded_primary_n == int(locked["a_target_case_reports"])
                and n_reports == int(locked["background_reports_n"])
            )
            reconciliation.append(
                {
                    "database": db,
                    "drug": drug,
                    "target_reports_rebuilt_n": loaded_target_n,
                    "target_reports_locked_n": int(locked["target_reports_n"]),
                    "primary_outcome_reports_rebuilt_n": loaded_primary_n,
                    "primary_outcome_reports_locked_n": int(locked["a_target_case_reports"]),
                    "background_reports_rebuilt_n": n_reports,
                    "background_reports_locked_n": int(locked["background_reports_n"]),
                    "reconciles_stage5": ok,
                }
            )
            if not ok:
                raise RuntimeError(f"Stage 5 reconciliation failed for {db}/{drug}")

            for outcome_key, pt in RENAL_PTS.items():
                row = s5.calc(
                    db,
                    drug,
                    outcome_key,
                    n_reports,
                    exposed,
                    events[outcome_key],
                    "exploratory_renal_preferred_term",
                )
                row["preferred_term"] = pt
                row[
                    "multiple_testing_family"
                ] = "All estimable Stage 5 secondary hypotheses plus renal-PT hypotheses; BY-adjusted"
                row[
                    "exact_95CI_interpretation"
                ] = "Nominal exact interval; BY q-value addresses exploratory multiplicity"
                row["_p"] = float(row["_p"])
                renal_rows.append(row)

                renal_cases = exposed & events[outcome_key]
                with_dka = renal_cases & dka
                with_dehydration = renal_cases & strict_dehydration
                with_dehydration_domain = renal_cases & dehydration_domain
                if not (
                    len(with_dka) <= len(renal_cases)
                    and len(with_dehydration) <= len(renal_cases)
                    and len(with_dehydration_domain) <= len(renal_cases)
                ):
                    raise RuntimeError(f"Invalid overlap arithmetic for {db}/{drug}/{pt}")
                renal_overlap_rows.append(
                    {
                        "database": db,
                        "drug": drug,
                        "renal_preferred_term": pt,
                        "target_reports_n": loaded_target_n,
                        "renal_PT_reports_n": len(renal_cases),
                        "renal_PT_and_DKA_reports_n": len(with_dka),
                        "renal_PT_and_Dehydration_PT_reports_n": len(with_dehydration),
                        "renal_PT_and_dehydration_electrolyte_domain_reports_n": len(
                            with_dehydration_domain
                        ),
                        "interpretation": "Report-level co-report counts; do not imply sequence, mechanism, or causation",
                    }
                )

            dka_cases = exposed & dka
            dka_dehydration = dka_cases & strict_dehydration
            dka_dehydration_domain = dka_cases & dehydration_domain
            if len(dka_dehydration) > len(dka_cases) or len(dka_dehydration_domain) > len(
                dka_cases
            ):
                raise RuntimeError(f"Invalid DKA/dehydration overlap arithmetic for {db}/{drug}")
            dka_dehydration_rows.append(
                {
                    "database": db,
                    "drug": drug,
                    "target_reports_n": loaded_target_n,
                    "DKA_reports_n": len(dka_cases),
                    "Dehydration_PT_reports_n": len(exposed & strict_dehydration),
                    "DKA_and_Dehydration_PT_reports_n": len(dka_dehydration),
                    "DKA_reports_with_Dehydration_PT_pct": round(
                        100 * len(dka_dehydration) / len(dka_cases), 2
                    )
                    if dka_cases
                    else "",
                    "dehydration_electrolyte_domain_reports_n": len(exposed & dehydration_domain),
                    "DKA_and_dehydration_electrolyte_domain_reports_n": len(dka_dehydration_domain),
                    "DKA_reports_with_dehydration_electrolyte_domain_pct": round(
                        100 * len(dka_dehydration_domain) / len(dka_cases), 2
                    )
                    if dka_cases
                    else "",
                    "interpretation": "Report-level co-report counts; do not imply event sequence, mechanism, or causation",
                }
            )

    if len(renal_rows) != 45 or len(dka_dehydration_rows) != 15 or len(renal_overlap_rows) != 45:
        raise RuntimeError(
            f"Unexpected output sizes: renal={len(renal_rows)}, DKA overlap={len(dka_dehydration_rows)}, renal overlap={len(renal_overlap_rows)}"
        )
    n_prior_secondary, n_combined_secondary = apply_by(
        renal_rows, "BY_adjusted_q_across_all_secondary_results"
    )
    renal_out = public_rows(renal_rows)
    write_csv(OUT / "阶段六_审稿人补充分析_renal_PT_ROR.csv", renal_out)
    write_csv(OUT / "阶段六_审稿人补充分析_DKA_Dehydration_overlap.csv", dka_dehydration_rows)
    write_csv(OUT / "阶段六_审稿人补充分析_renal_PT_co_report.csv", renal_overlap_rows)
    write_csv(OUT / "阶段六_审稿人补充分析_结果重建核验.csv", reconciliation)
    write_csv(OUT / "阶段六_审稿人补充分析_来源文件.csv", source_records)

    validation = {
        "analysis_date": datetime.now().isoformat(timespec="seconds"),
        "status": "completed_with_stage5_reconciliation",
        "design": "Exploratory report-level case-noncase ROR; database-specific; no pooling",
        "exposure_roles": {
            "FAERS": "primary suspect only",
            "Canada Vigilance": "native suspect role",
            "JADER": "suspected drug role",
        },
        "renal_preferred_terms": list(RENAL_PTS.values()),
        "multiplicity": "Benjamini-Yekutieli across 98 non-renal and 42 renal-PT estimable secondary comparisons (140 unique hypotheses total)",
        "nominal_confidence_intervals": "Exact conditional 95% CIs; not multiplicity-adjusted",
        "overlap_definition": "Report-level co-report; strict Dehydration PT and the broader dehydration/electrolyte domain are shown separately",
        "n_renal_comparisons": len(renal_rows),
        "n_estimable_renal_comparisons": sum(bool(r["estimable"]) for r in renal_rows),
        "n_prior_estimable_stage5_secondary_comparisons": n_prior_secondary,
        "n_combined_estimable_secondary_comparisons": n_combined_secondary,
        "n_dka_dehydration_rows": len(dka_dehydration_rows),
        "all_reconcile_stage5": all(r["reconciles_stage5"] for r in reconciliation),
        "all_overlap_checks_pass": True,
        "source_records": source_records,
        "source_fingerprint_scope": "Resolved paths, file sizes, and source modification times; raw-source SHA-256 was not computed",
    }
    (OUT / "阶段六_审稿人补充分析_核验.json").write_text(
        json.dumps(validation, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print(
        json.dumps(
            {
                "renal_comparisons": len(renal_rows),
                "renal_estimable": validation["n_estimable_renal_comparisons"],
                "DKA_dehydration_overlap_rows": len(dka_dehydration_rows),
                "stage5_reconciliations_passed": sum(
                    r["reconciles_stage5"] for r in reconciliation
                ),
                "outputs": [
                    "阶段六_审稿人补充分析_renal_PT_ROR.csv",
                    "阶段六_审稿人补充分析_DKA_Dehydration_overlap.csv",
                    "阶段六_审稿人补充分析_renal_PT_co_report.csv",
                    "阶段六_审稿人补充分析_结果重建核验.csv",
                    "阶段六_审稿人补充分析_来源文件.csv",
                    "阶段六_审稿人补充分析_核验.json",
                ],
            },
            ensure_ascii=False,
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
