"""Revision-stage post hoc JADER PT onset interval descriptions, aggregate output only."""
from __future__ import annotations
import argparse, csv, json, math, statistics, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
import audit_jader_clinical_availability as common


def linear_quantile(values, p):
    v = sorted(values)
    position = (len(v) - 1) * p
    lo = math.floor(position)
    hi = math.ceil(position)
    return v[lo] + (v[hi] - v[lo]) * (position - lo)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[3])
    p.add_argument("--source-dir", type=Path, required=True)
    p.add_argument("--stage5-script", type=Path, required=True)
    p.add_argument("--audit-dir", type=Path, default=Path(__file__).resolve().parent)
    p.add_argument("--out", type=Path, default=Path(__file__).resolve().parent)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    stage5 = a.stage5_script or a.root / "03_确认后补充分析/阶段五_主要分析/run_stage5_primary_analysis.py"
    s5 = common.module(stage5, "locked_jader_definitions")
    prior = json.loads((a.audit_dir / "data_audit.json").read_text(encoding="utf-8"))
    if prior["status"] != "audit_and_descriptions_completed" or not all(prior["checks"].values()):
        raise ValueError("Completed data audit required")
    rules = {
        "stage": "revision-stage post hoc exploratory description",
        "date_fixed": "2026-10-20",
        "snapshot_month_end": "2026-05-31",
        "strict_pair": "all corresponding target-suspect starts and same-PT onset dates complete valid, with one distinct date per field",
        "exclusions": [
            "incomplete, invalid or ambiguous start/onset dates",
            "negative interval",
            "any linked full start/end/onset date later than snapshot month end",
            "same drug-row end before start",
        ],
        "same_day_zero": "included",
        "onset_after_end": "flagged, not excluded",
        "quantiles": "linear interpolation at sorted index (n-1)*p; p=.25,.5,.75",
        "small_cells": "n<3 shown transparently; not generalized or compared",
        "output": "aggregate only; no IDs, exact dates or individual intervals",
    }
    (a.out / "TTO_rules_before_calculation.json").write_text(
        json.dumps(rules, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    meta = {}
    for source in prior["sources"]:
        if common.digest(a.source_dir / source["file"]) != source["sha256"]:
            raise ValueError("Source changed since audit")
    locked = {r["drug"]: r for r in common.read_csv(a.audit_dir / "cohort_reconciliation.csv")}
    audited = {
        (r["drug"], r["PT"]): r
        for r in common.read_csv(a.audit_dir / "PT_date_availability_by_drug.csv")
    }
    latest = {}
    for row in common.read_rows(a.source_dir / "demo202605.csv", "demo", meta):
        v = int(row[1])
        if row[0] not in latest or v > latest[row[0]]:
            latest[row[0]] = v
    reports = {s5.jkey(rid, v) for rid, v in latest.items()}
    target = {d: set() for d in s5.DRUG_TERMS}
    therapy = {d: defaultdict(list) for d in target}
    for row in common.read_rows(a.source_dir / "drug202605.csv", "drug", meta):
        key = s5.jkey(row[0], row[1])
        if key not in reports or row[3] != "被疑薬":
            continue
        d = s5.match_drug(row[4], row[5])
        if d:
            target[d].add(key)
            therapy[d][key].append((row[7], row[8]))
    onsets = {pt: defaultdict(list) for pt in common.PRIMARY}
    for row in common.read_rows(a.source_dir / "reac202605.csv", "reac", meta):
        key = s5.jkey(row[0], row[1])
        pt = s5.JP_TO_EN.get(row[3], "")
        if key in reports and pt in onsets:
            onsets[pt][key].append(row[5])
    rows = []
    checks = []
    for d, ids in target.items():
        primary = set.union(*(set(v) for v in onsets.values())) & ids
        if len(ids) != int(locked[d]["target_suspect_reports_n"]) or len(primary) != int(
            locked[d]["primary_reports_n"]
        ):
            raise ValueError("Frozen cohort mismatch")
        for pt in common.PRIMARY:
            ptids = ids & set(onsets[pt])
            values = []
            reason = defaultdict(int)
            after_end = 0
            after_latest = 0
            same_day = 0
            with_end = 0
            for key in ptids:
                drug_dates = therapy[d][key]
                ss = [common.date_status(x[0]) for x in drug_dates]
                ee = [common.date_status(x[1]) for x in drug_dates]
                oo = [common.date_status(v) for v in onsets[pt][key]]
                starts = {v for st, v in ss if st == "full_date_valid"}
                events = {v for st, v in oo if st == "full_date_valid"}
                ends = {v for st, v in ee if st == "full_date_valid"}
                strict = (
                    bool(ss)
                    and bool(oo)
                    and all(st == "full_date_valid" for st, v in ss + oo)
                    and len(starts) == len(events) == 1
                )
                if not strict:
                    reason["incomplete_invalid_or_ambiguous_dates_n"] += 1
                    continue
                start = next(iter(starts))
                event = next(iter(events))
                delta = (event - start).days
                if any(v and v > common.SNAPSHOT_MONTH_END for st, v in ss + ee + oo):
                    reason["postsnapshot_date_n"] += 1
                    continue
                if any(s and e and e < s for (_, s), (_, e) in zip(ss, ee)):
                    reason["end_before_start_n"] += 1
                    continue
                if delta < 0:
                    reason["negative_interval_n"] += 1
                    continue
                values.append(delta)
                same_day += delta == 0
                with_end += bool(ends)
                after_end += any(event > e for e in ends)
                after_latest += bool(ends) and event > max(ends)
            n = len(values)
            denom = len(ptids)
            reference = audited[(d, pt)]
            if n != int(
                reference["strict_pair_nonnegative_no_postsnapshot_date_no_end_conflict_reports_n"]
            ) or denom != int(reference["PT_reports_n"]):
                raise ValueError("TTO count does not match audited eligibility")
            r = {
                "drug": d,
                "PT": pt,
                "all_PT_reports_n": denom,
                "eligible_date_pairs_n": n,
                "excluded_or_unavailable_n": denom - n,
                "eligible_pct": 100 * n / denom if denom else "",
                "excluded_or_unavailable_pct": 100 * (denom - n) / denom if denom else "",
                "median_days": linear_quantile(values, 0.5) if n else "",
                "Q1_days_linear": linear_quantile(values, 0.25) if n else "",
                "Q3_days_linear": linear_quantile(values, 0.75) if n else "",
                "min_days": min(values) if n else "",
                "max_days": max(values) if n else "",
                "same_day_0_n": same_day,
                "eligible_with_any_full_end_date_n": with_end,
                "onset_after_any_full_end_date_n": after_end,
                "onset_after_latest_full_end_date_n": after_latest,
                **{
                    k: reason[k]
                    for k in [
                        "incomplete_invalid_or_ambiguous_dates_n",
                        "postsnapshot_date_n",
                        "end_before_start_n",
                        "negative_interval_n",
                    ]
                },
                "scope": "post hoc descriptive recorded interval; no comparison or causal inference",
            }
            rows.append(r)
            concordant = (
                (
                    all(
                        abs(
                            linear_quantile(values, p)
                            - float(np.quantile(values, p, method="linear"))
                        )
                        < 1e-12
                        for p in [0.25, 0.5, 0.75]
                    )
                    and linear_quantile(values, 0.5) == statistics.median(values)
                )
                if n
                else True
            )
            checks.append(
                {
                    "drug": d,
                    "PT": pt,
                    "audited_denominator_and_eligible_count_match": True,
                    "independent_quantiles_match": bool(concordant),
                    "partition_matches": n + sum(reason.values()) == denom,
                    "all_intervals_nonnegative": all(x >= 0 for x in values),
                }
            )
    if not all(all(v for k, v in r.items() if k not in {"drug", "PT"}) for r in checks):
        raise ValueError("Independent aggregate verification failed")
    common.save_csv(a.out / "JADER_PT_TTO_descriptive.csv", rows)
    common.save_csv(a.out / "JADER_PT_TTO_count_and_quantile_verification.csv", checks)
    manifest = {
        "status": "completed",
        "audit_required_and_passed": True,
        "rules_sha256": common.digest(a.out / "TTO_rules_before_calculation.json"),
        "source_audit_sha256": common.digest(a.audit_dir / "data_audit.json"),
        "date_availability_input_sha256": common.digest(
            a.audit_dir / "PT_date_availability_by_drug.csv"
        ),
        "stage5_definitions_sha256": common.digest(stage5),
        "source_files_sha256": {r["file"]: r["sha256"] for r in prior["sources"]},
        "python_version": sys.version.split()[0],
        "numpy_version": np.__version__,
        "verification": "exact cohort/eligibility count reconciliation; manual linear order-statistic interpolation vs NumPy linear and statistics.median, all in memory",
        "no_individual_identifiers_dates_or_intervals_exported": True,
        "no_effect_model_or_between_drug_test": True,
        "nonempty_cells_n": sum(r["all_PT_reports_n"] > 0 for r in rows),
        "all_cells_including_zero_denominators_n": len(rows),
    }
    (a.out / "JADER_PT_TTO_verification.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "status": "completed",
                "nonempty_results": [r for r in rows if r["all_PT_reports_n"]],
                "all_verification_passed": True,
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
