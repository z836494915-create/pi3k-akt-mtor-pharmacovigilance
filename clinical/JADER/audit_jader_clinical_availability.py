"""Read-only source audit; aggregate outputs only. No effect estimates or onset summaries."""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, json, re, unicodedata
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

AS_OF = date(2026, 10, 20)
SNAPSHOT_MONTH_END = date(2026, 5, 31)
PRIMARY = ["Diabetic ketoacidosis", "Ketoacidosis", "Metabolic acidosis"]
UNKNOWN = {"", "不明", "不明な病歴", "不明な原疾患", "unknown", "Unknown"}
DICTIONARY = {
    "recorded_diabetes_term": [
        "糖尿病",
        "1型糖尿病",
        "2型糖尿病",
        "インスリン依存性糖尿病",
        "インスリン非依存性糖尿病",
    ],
    "recorded_obesity_term": ["肥満", "肥満症", "病的肥満"],
}
DATE_STATUSES = [
    "blank",
    "unknown",
    "year_valid",
    "year_month_valid",
    "full_date_valid",
    "invalid_date_or_format",
]
HEADERS = {
    "demo": ["識別番号", "報告回数", "性別", "年齢", "体重", "身長", "報告年度・四半期", "状況", "報告の種類", "報告者の資格", "E2B"],
    "hist": ["識別番号", "報告回数", "原疾患等連番", "原疾患等"],
    "drug": [
        "識別番号",
        "報告回数",
        "医薬品連番",
        "医薬品の関与",
        "医薬品（一般名）",
        "医薬品（販売名）",
        "経路",
        "投与開始日",
        "投与終了日",
        "投与量",
        "投与単位",
        "分割投与回数",
        "使用理由",
        "医薬品の処置",
        "再投与による再発の有無",
        "リスク区分等",
    ],
    "reac": ["識別番号", "報告回数", "有害事象連番", "有害事象", "転帰", "有害事象の発現日"],
}


def read_rows(path, kind, audit):
    with path.open(encoding="cp932", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        if header != HEADERS[kind]:
            raise ValueError(f"{kind}: source schema mismatch")
        audit[kind + "_header"] = header
        for row in reader:
            if len(row) != len(header):
                raise ValueError(f"{kind}: unexpected column count")
            yield [x.strip() for x in row]


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def save_csv(path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def normalized(s):
    return unicodedata.normalize("NFKC", s.strip())


def value_status(s):
    if not s.strip():
        return "blank"
    if s.strip() in UNKNOWN:
        return "unknown"
    return "recorded_category"


def date_status(s):
    s = s.strip()
    if not s:
        return "blank", None
    if s in UNKNOWN:
        return "unknown", None
    if not re.fullmatch(r"\d{4}(?:\d{2})?(?:\d{2})?", s):
        return "invalid_date_or_format", None
    try:
        y = int(s[:4])
        if len(s) == 4:
            date(y, 1, 1)
            return "year_valid", None
        m = int(s[4:6])
        if len(s) == 6:
            date(y, m, 1)
            return "year_month_valid", None
        return "full_date_valid", date(y, m, int(s[6:8]))
    except ValueError:
        return "invalid_date_or_format", None


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[3])
    p.add_argument("--source", "--source-dir", dest="source", type=Path, required=True)
    p.add_argument("--stage5-script", type=Path, required=True)
    p.add_argument("--cohort-reference", type=Path, required=True)
    p.add_argument("--out", type=Path, default=Path(__file__).resolve().parent)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    s5_path = a.stage5_script or a.root / "03_确认后补充分析/阶段五_主要分析/run_stage5_primary_analysis.py"
    ref_path = (
        a.cohort_reference
        or a.root / "03_确认后补充分析/阶段八_核验后修复/External/external_stage5_reconciliation.csv"
    )
    s5 = module(s5_path, "stage5_definitions")
    locked = {r["drug"]: r for r in read_csv(ref_path) if r["database"] == "JADER"}
    audit = {
        "as_of": str(AS_OF),
        "design": "source data audit and descriptive availability only",
        "sources": [],
        "stage5_definition_file_sha256": digest(s5_path),
        "cohort_reference_sha256": digest(ref_path),
        "linkage": "highest numeric DEMO version per ID; exact ID and report version on all tables",
        "drug_role": "被疑薬",
        "clinical_history_interpretation": "recorded underlying-condition terms, not validated pretreatment comorbidity",
        "clinical_validation": False,
        "effect_models_run": False,
        "TTO_distribution_computed": False,
    }
    rules = {
        "as_of": str(AS_OF),
        "snapshot_month_end": str(SNAPSHOT_MONTH_END),
        "primary_PT": PRIMARY,
        "drug_terms": s5.DRUG_TERMS,
        "history_dictionary": DICTIONARY,
        "date_statuses": DATE_STATUSES,
        "strict_unique_pair": "all target-suspect start rows and all matching PT onset rows full valid; exactly one distinct start and onset date",
        "partial_date_imputation": False,
        "height_weight": "original categories only; no BMI derivation",
        "outputs": "aggregate counts only; no report IDs, exact dates, individual intervals",
    }
    (a.out / "rules_before_audit.json").write_text(
        json.dumps(rules, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    for kind in HEADERS:
        f = a.source / f"{kind}202605.csv"
        audit["sources"].append({"file": f.name, "bytes": f.stat().st_size, "sha256": digest(f)})
    latest = {}
    seen = set()
    dups = 0
    demo_n = 0
    for row in read_rows(a.source / "demo202605.csv", "demo", audit):
        demo_n += 1
        key = s5.jkey(row[0], row[1])
        dups += key in seen
        seen.add(key)
        v = int(row[1])
        if row[0] not in latest or v > latest[row[0]]:
            latest[row[0]] = v
    if dups:
        raise ValueError("DEMO duplicate ID-version rows require review")
    reports = {s5.jkey(r, v) for r, v in latest.items()}
    del seen
    target = {d: set() for d in s5.DRUG_TERMS}
    starts = {d: defaultdict(list) for d in target}
    ends = {d: defaultdict(list) for d in target}
    sequences = {d: defaultdict(set) for d in target}
    role_counts = Counter()
    drug_n = 0
    drug_eligible = 0
    drug_duplicate_seq = 0
    drug_ambiguous_seq = 0
    drug_seen = {}
    for row in read_rows(a.source / "drug202605.csv", "drug", audit):
        drug_n += 1
        role_counts[row[3]] += 1
        key = s5.jkey(row[0], row[1])
        if key not in reports or row[3] != "被疑薬":
            continue
        d = s5.match_drug(row[4], row[5])
        if not d:
            continue
        drug_eligible += 1
        target[d].add(key)
        starts[d][key].append(row[7])
        ends[d][key].append(row[8])
        sequences[d][key].add(row[2])
        sk = (key, row[2])
        sig = tuple(row)
        if sk in drug_seen:
            drug_duplicate_seq += 1
            drug_ambiguous_seq += drug_seen[sk] != sig
        else:
            drug_seen[sk] = sig
    primary_pts = {pt: set() for pt in PRIMARY}
    onsets = {pt: defaultdict(list) for pt in PRIMARY}
    reac_n = 0
    orphan_reac = 0
    reaction_duplicate_seq = 0
    reaction_ambiguous_seq = 0
    reaction_seen = {}
    for row in read_rows(a.source / "reac202605.csv", "reac", audit):
        reac_n += 1
        key = s5.jkey(row[0], row[1])
        if key not in reports:
            orphan_reac += 1
            continue
        en = s5.JP_TO_EN.get(row[3], "")
        if en not in PRIMARY:
            continue
        primary_pts[en].add(key)
        onsets[en][key].append(row[5])
        sk = (key, row[2])
        sig = tuple(row)
        if sk in reaction_seen:
            reaction_duplicate_seq += 1
            reaction_ambiguous_seq += reaction_seen[sk] != sig
        else:
            reaction_seen[sk] = sig
    primary = set.union(*primary_pts.values())
    cohorts = {d: ids & primary for d, ids in target.items()}
    selected = set.union(*cohorts.values())
    reconciliation = []
    for d in target:
        old = locked[d]
        passed = (
            len(reports) == int(old["background_reports_n"])
            and len(target[d]) == int(old["target_reports_n"])
            and len(cohorts[d]) == int(old["primary_reports_n"])
        )
        reconciliation.append(
            {
                "drug": d,
                "all_database_reports_n": len(reports),
                "target_suspect_reports_n": len(target[d]),
                "primary_reports_n": len(cohorts[d]),
                "locked_counts_match": passed,
            }
        )
        if not passed:
            raise ValueError(f"Locked count mismatch: {d}")
    save_csv(a.out / "cohort_reconciliation.csv", reconciliation)
    print("Locked JADER counts reconciled; auditing available fields.", flush=True)
    hist_terms = defaultdict(list)
    hist_total = 0
    hist_versions = set()
    hist_ids = set()
    hist_seen = {}
    hist_duplicate_seq = 0
    hist_exact_dups = 0
    hist_conflicts = 0
    hist_orphan = 0
    hist_old = 0
    hist_blank_terms = 0
    hist_row_lengths = Counter()
    for row in read_rows(a.source / "hist202605.csv", "hist", audit):
        hist_total += 1
        key = s5.jkey(row[0], row[1])
        hist_versions.add(key)
        hist_ids.add(row[0])
        hist_blank_terms += not bool(row[3])
        sk = (key, row[2])
        term = row[3]
        if sk in hist_seen:
            hist_duplicate_seq += 1
            hist_exact_dups += term == hist_seen[sk]
            hist_conflicts += term != hist_seen[sk]
        else:
            hist_seen[sk] = term
        if row[0] not in latest:
            hist_orphan += 1
        elif int(row[1]) != latest[row[0]]:
            hist_old += 1
        if key in selected:
            hist_terms[key].append(term)
    del hist_seen, drug_seen, reaction_seen
    demo = {}
    for row in read_rows(a.source / "demo202605.csv", "demo", audit):
        key = s5.jkey(row[0], row[1])
        if key in selected:
            demo[key] = row
    audit.update(
        {
            "demo_rows_n": demo_n,
            "latest_reports_n": len(reports),
            "superseded_demo_rows_n": demo_n - len(reports),
            "demo_duplicate_ID_version_rows_n": dups,
            "target_suspect_drug_rows_n": drug_eligible,
            "target_suspect_duplicate_sequence_rows_n": drug_duplicate_seq,
            "target_suspect_conflicting_sequence_rows_n": drug_ambiguous_seq,
            "drug_rows_n": drug_n,
            "drug_roles": dict(role_counts),
            "reaction_rows_n": reac_n,
            "reaction_nonretained_ID_version_rows_n": orphan_reac,
            "primary_reaction_duplicate_sequence_rows_n": reaction_duplicate_seq,
            "primary_reaction_conflicting_sequence_rows_n": reaction_ambiguous_seq,
            "hist_rows_n": hist_total,
            "hist_distinct_IDs_n": len(hist_ids),
            "hist_distinct_ID_version_n": len(hist_versions),
            "hist_duplicate_ID_version_sequence_rows_n": hist_duplicate_seq,
            "hist_exact_duplicate_rows_n": hist_exact_dups,
            "hist_conflicting_same_sequence_rows_n": hist_conflicts,
            "hist_orphan_ID_rows_n": hist_orphan,
            "hist_nonlatest_version_rows_n": hist_old,
            "hist_blank_term_rows_n": hist_blank_terms,
            "target_primary_unique_reports_n": len(selected),
            "source_counts_match_locked_results": True,
        }
    )
    dictionary_rows = [
        {
            "category": cat,
            "exact_normalized_term": term,
            "interpretation": "recorded term only; not confirmed pretreatment condition",
        }
        for cat, terms in DICTIONARY.items()
        for term in terms
    ]
    save_csv(a.out / "history_dictionary_fixed_before_counts.csv", dictionary_rows)
    hist_summary = []
    term_rows = []
    anthro_rows = []
    patterns = []
    date_reports = []
    date_precision = []
    date_patterns = []
    for d, ids in cohorts.items():
        terms_count = Counter()
        hs = Counter()
        cat_ids = {cat: set() for cat in DICTIONARY}
        pc = Counter()
        ac = Counter()
        for key in ids:
            terms = hist_terms.get(key, [])
            known = {t for t in terms if t not in UNKNOWN}
            hs["hist_any_row_reports_n"] += bool(terms)
            hs["hist_known_nonempty_term_reports_n"] += bool(known)
            hs["hist_only_blank_unknown_reports_n"] += bool(terms) and not known
            hs["hist_no_row_reports_n"] += not terms
            hs["hist_raw_rows_n"] += len(terms)
            hs["hist_unique_term_occurrences_n"] += len(set(terms))
            for t in set(terms):
                terms_count[t] += 1
            for cat, dictionary in DICTIONARY.items():
                if any(normalized(t) in dictionary for t in known):
                    cat_ids[cat].add(key)
            ws = value_status(demo[key][4])
            hh = value_status(demo[key][5])
            ac[("weight", ws, demo[key][4])] += 1
            ac[("height", hh, demo[key][5])] += 1
            start_ok = all(date_status(s)[0] == "full_date_valid" for s in starts[d][key])
            anyonset = any(
                any(date_status(v)[0] == "full_date_valid" for v in onsets[pt].get(key, []))
                for pt in PRIMARY
            )
            pc[
                (
                    not bool(known),
                    ws != "recorded_category",
                    hh != "recorded_category",
                    not start_ok,
                    not anyonset,
                )
            ] += 1
        hist_summary.append(
            {
                "drug": d,
                "primary_reports_n": len(ids),
                **{
                    k: hs[k]
                    for k in [
                        "hist_any_row_reports_n",
                        "hist_known_nonempty_term_reports_n",
                        "hist_only_blank_unknown_reports_n",
                        "hist_no_row_reports_n",
                        "hist_raw_rows_n",
                        "hist_unique_term_occurrences_n",
                    ]
                },
                **{cat + "_reports_n": len(v) for cat, v in cat_ids.items()},
                "hist_known_term_missing_or_unknown_n": len(ids)
                - hs["hist_known_nonempty_term_reports_n"],
                "hist_known_term_missing_or_unknown_pct": 100
                * (len(ids) - hs["hist_known_nonempty_term_reports_n"])
                / len(ids)
                if ids
                else "",
            }
        )
        for term, n in sorted(terms_count.items()):
            cats = [cat for cat, dd in DICTIONARY.items() if normalized(term) in dd]
            term_rows.append(
                {
                    "drug": d,
                    "primary_reports_n": len(ids),
                    "original_hist_term": term,
                    "normalized_term": normalized(term),
                    "reports_with_term_n": n,
                    "dictionary_category": ";".join(cats),
                    "term_status": "blank_or_unknown" if term in UNKNOWN else "recorded_term",
                    "interpretation": "original recorded term; no verified temporal or clinical interpretation",
                }
            )
        for (field, status, value), n in sorted(ac.items()):
            anthro_rows.append(
                {
                    "drug": d,
                    "primary_reports_n": len(ids),
                    "field": field,
                    "category_status": status,
                    "original_category": value,
                    "reports_n": n,
                }
            )
        for combo, n in sorted(pc.items()):
            patterns.append(
                {
                    "drug": d,
                    "primary_reports_n": len(ids),
                    "history_known_term_unavailable": combo[0],
                    "weight_unavailable": combo[1],
                    "height_unavailable": combo[2],
                    "all_target_start_rows_full_valid_unavailable": combo[3],
                    "any_primary_onset_full_valid_unavailable": combo[4],
                    "reports_n": n,
                }
            )
        for pt in PRIMARY:
            ptids = ids & primary_pts[pt]
            summary = Counter()
            prec = Counter()
            pats = Counter()
            for key in ptids:
                sv = starts[d][key]
                ev = ends[d][key]
                ov = onsets[pt][key]
                sf = [date_status(x) for x in sv]
                ef = [date_status(x) for x in ev]
                of = [date_status(x) for x in ov]
                ss = {v for status, v in sf if status == "full_date_valid"}
                es = {v for status, v in ef if status == "full_date_valid"}
                os = {v for status, v in of if status == "full_date_valid"}
                for field, values in [
                    ("target_drug_start", sf),
                    ("target_drug_end", ef),
                    ("PT_onset", of),
                ]:
                    for status, val in values:
                        prec[(field, status)] += 1
                    summary[field + "_any_full_valid_reports_n"] += any(
                        status == "full_date_valid" for status, val in values
                    )
                    summary[field + "_all_rows_full_valid_reports_n"] += bool(values) and all(
                        status == "full_date_valid" for status, val in values
                    )
                    summary[field + "_any_future_full_date_reports_n"] += any(
                        val and val > AS_OF for status, val in values
                    )
                    summary[field + "_any_after_snapshot_month_end_reports_n"] += any(
                        val and val > SNAPSHOT_MONTH_END for status, val in values
                    )
                summary["target_drug_multiple_rows_reports_n"] += len(sv) > 1
                summary["target_drug_multiple_sequences_reports_n"] += len(sequences[d][key]) > 1
                summary["target_start_multiple_distinct_recorded_values_reports_n"] += (
                    len({x for x in sv if x}) > 1
                )
                summary["PT_onset_multiple_rows_reports_n"] += len(ov) > 1
                summary["PT_onset_multiple_distinct_recorded_values_reports_n"] += (
                    len({x for x in ov if x}) > 1
                )
                summary["target_start_multiple_distinct_full_dates_reports_n"] += len(ss) > 1
                summary["PT_onset_multiple_distinct_full_dates_reports_n"] += len(os) > 1
                end_conflict = any(sv1 and ev1 and ev1 < sv1 for (_, sv1), (_, ev1) in zip(sf, ef))
                summary["any_target_end_before_start_reports_n"] += end_conflict
                summary["any_full_date_combination_negative_interval_reports_n"] += any(
                    o < s for s in ss for o in os
                )
                summary["any_full_date_onset_after_end_reports_n"] += any(
                    o > e for o in os for e in es
                )
                strict = (
                    bool(sf)
                    and bool(of)
                    and all(t == "full_date_valid" for t, v in sf + of)
                    and len(ss) == len(os) == 1
                )
                summary["strict_unique_full_date_pair_reports_n"] += strict
                if strict:
                    sd = next(iter(ss))
                    od = next(iter(os))
                    delta = (od - sd).days
                    status = "negative" if delta < 0 else "same_day" if delta == 0 else "positive"
                    summary["strict_pair_" + status + "_reports_n"] += 1
                    future = sd > AS_OF or od > AS_OF
                    summary["strict_pair_future_date_reports_n"] += future
                    summary["strict_pair_after_snapshot_month_end_reports_n"] += (
                        sd > SNAPSHOT_MONTH_END or od > SNAPSHOT_MONTH_END
                    )
                    summary["strict_pair_end_before_start_reports_n"] += end_conflict
                    summary["strict_pair_onset_after_any_full_end_reports_n"] += any(
                        od > e for e in es
                    )
                    summary["strict_pair_onset_after_latest_full_end_reports_n"] += bool(
                        es
                    ) and od > max(es)
                    summary["strict_pair_nonnegative_no_future_no_end_conflict_reports_n"] += (
                        delta >= 0 and not future and not end_conflict
                    )
                    summary[
                        "strict_pair_nonnegative_no_postsnapshot_date_no_end_conflict_reports_n"
                    ] += (
                        delta >= 0
                        and sd <= SNAPSHOT_MONTH_END
                        and od <= SNAPSHOT_MONTH_END
                        and all(e <= SNAPSHOT_MONTH_END for e in es)
                        and not end_conflict
                    )
                pats[(bool(ss), bool(os), len(ss) > 1, len(os) > 1, strict)] += 1
            keys = [
                "target_drug_multiple_rows_reports_n",
                "target_drug_multiple_sequences_reports_n",
                "target_start_multiple_distinct_recorded_values_reports_n",
                "PT_onset_multiple_rows_reports_n",
                "PT_onset_multiple_distinct_recorded_values_reports_n",
                "target_start_multiple_distinct_full_dates_reports_n",
                "PT_onset_multiple_distinct_full_dates_reports_n",
                "any_target_end_before_start_reports_n",
                "any_full_date_combination_negative_interval_reports_n",
                "any_full_date_onset_after_end_reports_n",
                "strict_unique_full_date_pair_reports_n",
                "strict_pair_negative_reports_n",
                "strict_pair_same_day_reports_n",
                "strict_pair_positive_reports_n",
                "strict_pair_future_date_reports_n",
                "strict_pair_after_snapshot_month_end_reports_n",
                "strict_pair_end_before_start_reports_n",
                "strict_pair_onset_after_any_full_end_reports_n",
                "strict_pair_onset_after_latest_full_end_reports_n",
                "strict_pair_nonnegative_no_future_no_end_conflict_reports_n",
                "strict_pair_nonnegative_no_postsnapshot_date_no_end_conflict_reports_n",
            ]
            keys += [
                field + suffix
                for field in ["target_drug_start", "target_drug_end", "PT_onset"]
                for suffix in [
                    "_any_full_valid_reports_n",
                    "_all_rows_full_valid_reports_n",
                    "_any_future_full_date_reports_n",
                    "_any_after_snapshot_month_end_reports_n",
                ]
            ]
            date_reports.append(
                {"drug": d, "PT": pt, "PT_reports_n": len(ptids), **{k: summary[k] for k in keys}}
            )
            for field in ["target_drug_start", "target_drug_end", "PT_onset"]:
                for status in DATE_STATUSES:
                    date_precision.append(
                        {
                            "drug": d,
                            "PT": pt,
                            "PT_reports_n": len(ptids),
                            "field": field,
                            "date_status": status,
                            "rows_n": prec[(field, status)],
                        }
                    )
            for combo, n in sorted(pats.items()):
                date_patterns.append(
                    {
                        "drug": d,
                        "PT": pt,
                        "start_any_full_valid": combo[0],
                        "onset_any_full_valid": combo[1],
                        "start_multiple_full_dates": combo[2],
                        "onset_multiple_full_dates": combo[3],
                        "strict_unique_pair": combo[4],
                        "reports_n": n,
                    }
                )
    save_csv(a.out / "history_availability_by_drug.csv", hist_summary)
    save_csv(a.out / "history_original_terms_primary_reports.csv", term_rows)
    save_csv(a.out / "height_weight_original_categories.csv", anthro_rows)
    save_csv(a.out / "clinical_availability_missing_patterns.csv", patterns)
    save_csv(a.out / "PT_date_availability_by_drug.csv", date_reports)
    save_csv(a.out / "PT_date_precision_row_counts.csv", date_precision)
    save_csv(a.out / "PT_date_missing_patterns.csv", date_patterns)
    checks = {
        "all_primary_counts_reconcile": all(r["locked_counts_match"] for r in reconciliation),
        "history_categories_sum_to_denominator": all(
            r["hist_known_nonempty_term_reports_n"]
            + r["hist_only_blank_unknown_reports_n"]
            + r["hist_no_row_reports_n"]
            == r["primary_reports_n"]
            for r in hist_summary
        ),
        "strict_pair_categories_sum": all(
            r["strict_pair_negative_reports_n"]
            + r["strict_pair_same_day_reports_n"]
            + r["strict_pair_positive_reports_n"]
            == r["strict_unique_full_date_pair_reports_n"]
            for r in date_reports
        ),
        "clinical_pattern_counts_reconcile": all(
            sum(r["reports_n"] for r in patterns if r["drug"] == d) == len(ids)
            for d, ids in cohorts.items()
        ),
        "PT_date_pattern_counts_reconcile": all(
            sum(
                r["reports_n"]
                for r in date_patterns
                if r["drug"] == row["drug"] and r["PT"] == row["PT"]
            )
            == row["PT_reports_n"]
            for row in date_reports
        ),
    }
    if not all(checks.values()):
        raise ValueError("Aggregate reconciliation failed")
    audit["checks"] = checks
    audit["status"] = "audit_and_descriptions_completed"
    audit["rules_sha256"] = digest(a.out / "rules_before_audit.json")
    (a.out / "data_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    dd = {
        "unit": "unique retained report, not validated unique patient; a report can contribute to more than one drug/PT",
        "history_absence": "No HIST row / blank or unknown term is unrecorded information, not absence of disease",
        "history_terms": "terms counted once per report; multiple terms can coexist; dictionary is exact normalized mapping only",
        "date_rows": "DRUG dates counted within each drug/PT subset; same drug row can be counted in several PT subsets; do not sum as patients",
        "strict_unique_full_date_pair_reports_n": rules["strict_unique_pair"],
        "negative_interval": "onset before treatment start among available date combinations; descriptive quality flag only",
        "same_day": "day-level equality does not establish within-day exposure-before-event",
        "future_date": "full date later than 2026-10-20",
        "after_snapshot_month_end": "full date later than 2026-05-31, a separate snapshot-time inconsistency flag",
        "onset_after_end": "reported descriptively; not excluded because reactions may occur after discontinuation",
        "height_weight": "original categories only; no units or category midpoints inferred; not treatment-baseline measures",
        "zero_denominator": "zero-report subsets shown as zero counts, missingness percentage left empty",
    }
    (a.out / "data_dictionary.json").write_text(
        json.dumps(dd, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "status": audit["status"],
                "history_summary": hist_summary,
                "date_summary_nonempty": [r for r in date_reports if r["PT_reports_n"]],
                "checks": checks,
            },
            ensure_ascii=False,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
