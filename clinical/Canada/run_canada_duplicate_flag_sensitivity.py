"""Exploratory exclusion of officially Duplicate-associated Canada reports.

No Linked records are excluded. This is report-set sensitivity, not patient-level
deduplication. The original Stage 5 analysis and its multiplicity families remain
unchanged. Five new primary comparisons constitute a separate exploratory BY family.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib.util
import json
import math
from collections import Counter
from datetime import datetime
from pathlib import Path


def write_csv(path, rows):
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def rows(path, ncols):
    with path.open(encoding="utf-8-sig", newline="") as f:
        for r in csv.reader(f, delimiter="$"):
            if not r:
                continue
            if len(r) != ncols:
                raise ValueError("Source schema changed")
            yield [x.strip() for x in r]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[3])
    p.add_argument("--stage5", type=Path, required=True)
    p.add_argument("--canada", type=Path, required=True)
    p.add_argument("--out", type=Path, default=Path(__file__).resolve().parent)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    s5dir = a.stage5 or a.root / "03_确认后补充分析" / "阶段五_主要分析"
    spec = importlib.util.spec_from_file_location(
        "locked_stage5", s5dir / "run_stage5_primary_analysis.py"
    )
    s5 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(s5)
    with (s5dir / "阶段五_主要结局ROR_Fisher_Holm.csv").open(encoding="utf-8-sig", newline="") as f:
        locked = {r["drug"]: r for r in csv.DictReader(f) if r["database"] == "Canada Vigilance"}
    plan = {
        "date_fixed_before_calculation": datetime.now().isoformat(timespec="seconds"),
        "analysis": "Revision-stage exploratory Duplicate-flag exclusion sensitivity",
        "trigger": "Evaluate the effect of excluding officially Duplicate-associated reports; connected components are not adjudicated patients.",
        "population": "The same complete 2026-02-28 Canada extract as locked Stage 5.",
        "exposure": "The original five drugs under native Suspect matching.",
        "outcome": s5.OUTCOMES["primary_union"],
        "exclusion": "For every report_links row with RECORD_TYPE_ENG=Duplicate, exclude its source REPORT_ID if retained, including sources whose linked number is missing/unmatched. Also exclude each retained destination reached by an exact REPORT_LINK_NO=REPORT_NO mapping. Apply the single exclusion set to the entire report population, drug groups and outcome groups.",
        "not_excluded": "Linked rows do not generate exclusions; no inferred fuzzy link or undocumented identifier conversion is used.",
        "estimator": "Original cross-product ROR, conditional exact 95% CI and two-sided Fisher exact P.",
        "multiplicity": "The five new drug-primary comparisons form a separate exploratory Benjamini-Yekutieli family; retain nonestimable rows but adjust only finite estimable tests. Original Holm/BY families remain unchanged.",
        "interpretation": "A change in the reporting background after selective removal, not adjudicated patient deduplication, independent validation, incidence or causal risk.",
        "limitations": [
            "Flagged-report exclusion may alter clinical and reporting composition.",
            "Unflagged and cross-database duplicate reports may persist.",
            "Missing linked destinations cannot be recovered from the current extract.",
            "No representative is chosen from a connected component.",
        ],
        "schema_url": "https://www.canada.ca/en/health-canada/services/drugs-health-products/medeffect-canada/adverse-reaction-database/canada-vigilance-adverse-reaction-online-database-data-structures.html",
        "glossary_url": "https://www.canada.ca/en/health-canada/services/drugs-health-products/medeffect-canada/adverse-reaction-database/glossary.html",
        "protocol_timing": "Symmetric flag-exclusion protocol fixed after the source audit and before calculation.",
    }
    (a.out / "duplicate_sensitivity_plan_before_run.json").write_text(
        json.dumps(plan, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    report_set = set()
    no_to_id = {}
    for r in rows(a.canada / "reports.txt", 42):
        if r[0] in report_set or r[1] in no_to_id:
            raise ValueError("Report IDs/numbers no longer unique")
        report_set.add(r[0])
        no_to_id[r[1]] = r[0]
    print("Rebuilding locked target/outcome sets", flush=True)
    target = {d: set() for d in s5.DRUG_TERMS}
    for r in rows(a.canada / "report_drug.txt", 22):
        if r[1] in report_set and r[4].upper() == "SUSPECT":
            d = s5.match_drug(r[3])
            if d:
                target[d].add(r[1])
    pts = {x.lower() for x in s5.OUTCOMES["primary_union"]}
    cases = set()
    for r in rows(a.canada / "reactions.txt", 10):
        if r[1] in report_set and r[5].lower() in pts:
            cases.add(r[1])
    excluded_source = set()
    excluded_destination = set()
    link_counts = Counter()
    for r in rows(a.canada / "report_links.txt", 5):
        if r[2].lower() != "duplicate":
            continue
        link_counts["duplicate_rows"] += 1
        if r[1] in report_set:
            excluded_source.add(r[1])
        else:
            link_counts["source_not_retained"] += 1
        destination = no_to_id.get(r[4])
        if destination:
            excluded_destination.add(destination)
            link_counts["destination_exact_REPORT_NO"] += 1
            # Independent consistency: a REPORT_ID exact match must never refer
            # to a different report when an exact REPORT_NO mapping exists.
            if r[4] in report_set and r[4] != destination:
                raise ValueError("Ambiguous report-number versus report-ID mapping")
        elif r[4]:
            link_counts["destination_not_retained"] += 1
            if r[4] in report_set:
                raise ValueError(
                    "REPORT_ID-only destination observed: inspect schema before analysis"
                )
        else:
            link_counts["destination_empty"] += 1
    excluded = excluded_source | excluded_destination
    retained = report_set - excluded
    cases_new = cases & retained
    result = []
    checks = []
    cells = [
        "a_target_case_reports",
        "b_target_noncase_reports",
        "c_background_case_reports",
        "d_background_noncase_reports",
    ]
    for drug, exposed in target.items():
        old = locked[drug]
        rebuilt = s5.calc(
            "Canada Vigilance", drug, "primary_union", len(report_set), exposed, cases, "primary"
        )
        count_ok = all(
            int(old[k]) == rebuilt[k]
            for k in cells + ["target_reports_n", "background_reports_n", "total_case_reports_n"]
        )
        stats_ok = all(
            old[k] == rebuilt[k]
            for k in [
                "ROR_cross_product",
                "conditional_OR_MLE",
                "exact_conditional_95CI_lower",
                "exact_conditional_95CI_upper",
                "fisher_two_sided_p",
            ]
        )
        if not count_ok or not stats_ok:
            raise ValueError(f"Original locked analysis does not reconcile for {drug}")
        exposed_new = exposed & retained
        new = s5.calc(
            "Canada Vigilance",
            drug,
            "primary_union",
            len(retained),
            exposed_new,
            cases_new,
            "exploratory_duplicate_flag_exclusion",
        )
        removed_sets = [
            excluded & exposed & cases,
            (excluded & exposed) - cases,
            (excluded & cases) - exposed,
            excluded - (exposed | cases),
        ]
        for key, removed in zip(cells, removed_sets):
            if rebuilt[key] - new[key] != len(removed):
                raise ValueError("Independent 2x2 exclusion count mismatch")
            new["removed_" + key] = len(removed)
        if sum(len(s) for s in removed_sets) != len(excluded):
            raise ValueError("Removed cells do not sum to exclusion set")
        for k, v in old.items():
            if k not in [
                "database",
                "drug",
                "outcome",
                "outcome_class",
                "prespecified_primary_status",
            ]:
                new["original_" + k] = v
        new["removed_total_reports_n"] = len(excluded)
        new["removed_target_reports_n"] = len(exposed & excluded)
        new["removed_all_primary_reports_n"] = len(cases & excluded)
        new["removed_target_primary_reports_n"] = len(exposed & cases & excluded)
        new["exploratory_family"] = "five_canada_primary_Duplicate_flag_exclusion"
        new["is_patient_level_deduplication"] = False
        new["linked_rows_used_for_exclusion"] = False
        result.append(new)
        checks.append(
            {
                "drug": drug,
                "original_counts_match": count_ok,
                "original_statistics_match": stats_ok,
                "independent_removed_cell_counts_match": True,
            }
        )
    s5.adjust(result, "fdr_by", "exploratory_BY_q")
    # Independent BY implementation from the unrounded Fisher P values.
    est = [r for r in result if r["estimable"] and math.isfinite(r["_p"])]
    m = len(est)
    ordered = sorted(est, key=lambda r: r["_p"])
    harmonic = sum(1 / i for i in range(1, m + 1))
    running = 1.0
    independent = {}
    for rank in range(m, 0, -1):
        r = ordered[rank - 1]
        running = min(running, r["_p"] * m * harmonic / rank)
        independent[r["drug"]] = s5.fmt(running)
    for r in result:
        if r["estimable"]:
            if r["exploratory_BY_q"] != independent[r["drug"]]:
                raise ValueError("Independent BY calculation mismatch")
        else:
            r["exploratory_BY_q"] = ""
    write_csv(a.out / "duplicate_flag_exclusion_primary_results.csv", s5.clean(result))
    write_csv(a.out / "duplicate_flag_exclusion_reconciliation.csv", checks)
    validation = {
        "status": "completed_after_source_audit_and_fixed_protocol",
        "original_reports_n": len(report_set),
        "original_primary_reports_n": len(cases),
        "flagged_source_reports_n": len(excluded_source),
        "exact_linked_destination_reports_n": len(excluded_destination),
        "source_destination_overlap_n": len(excluded_source & excluded_destination),
        "excluded_union_reports_n": len(excluded),
        "excluded_primary_reports_n": len(excluded & cases),
        "retained_reports_n": len(retained),
        "retained_primary_reports_n": len(cases_new),
        "link_counts": dict(link_counts),
        "planned_contrasts": 5,
        "estimable_contrasts": m,
        "original_locked_counts_and_statistics_reproduced": True,
        "independent_2x2_cell_removal_check": True,
        "independent_BY_check": True,
        "original_analysis_outputs_modified": False,
        "patient_level_deduplication_claim": False,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "stage5_dependency_sha256": hashlib.sha256(
            (s5dir / "run_stage5_primary_analysis.py").read_bytes()
        ).hexdigest(),
    }
    (a.out / "duplicate_flag_exclusion_validation.json").write_text(
        json.dumps(validation, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(
        json.dumps({"validation": validation, "results": s5.clean(result)}, ensure_ascii=False),
        flush=True,
    )


if __name__ == "__main__":
    main()
