"""Read-only Canada field/duplicate-link audit; aggregates only, no new models.

Definitions fixed before running: original Stage 5 native-Suspect drug matching,
original three-PT union, one report per existing REPORT_ID/REPORT_NO. Numeric
fields are assessed as recorded, finite/positive, and accompanied by a recognized
unit; no BMI, baseline confounder, duration-to-onset or patient-level inference.
No records are removed based on duplicate links in this audit.
"""
from __future__ import annotations
import argparse
import csv
import gc
import hashlib
import importlib.util
import json
import math
import re
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

SCHEMA_URL = "https://www.canada.ca/en/health-canada/services/drugs-health-products/medeffect-canada/adverse-reaction-database/canada-vigilance-adverse-reaction-online-database-data-structures.html"
DURATION_UNITS = {
    "seconds",
    "second",
    "minutes",
    "minute",
    "hours",
    "hour",
    "days",
    "day",
    "weeks",
    "week",
    "months",
    "month",
    "years",
    "year",
}
WEIGHT_UNITS = {"kg", "kilogram", "kilograms", "g", "gram", "grams", "lb", "lbs", "pound", "pounds"}
HEIGHT_UNITS = {
    "cm",
    "centimeter",
    "centimeters",
    "centimetre",
    "centimetres",
    "m",
    "meter",
    "meters",
    "metre",
    "metres",
    "in",
    "inch",
    "inches",
    "ft",
    "foot",
    "feet",
}
WEIGHT_TO_KG = {
    **{k: 1 for k in ["kg", "kilogram", "kilograms"]},
    **{k: 0.001 for k in ["g", "gram", "grams"]},
    **{k: 0.45359237 for k in ["lb", "lbs", "pound", "pounds"]},
}
HEIGHT_TO_CM = {
    **{k: 1 for k in ["cm", "centimeter", "centimeters", "centimetre", "centimetres"]},
    **{k: 100 for k in ["m", "meter", "meters", "metre", "metres"]},
    **{k: 2.54 for k in ["in", "inch", "inches"]},
    **{k: 30.48 for k in ["ft", "foot", "feet"]},
}


def save_csv(path, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def read_rows(path, columns):
    with path.open(encoding="utf-8-sig", newline="") as f:
        for row in csv.reader(f, delimiter="$"):
            if not row:
                continue
            if len(row) != columns:
                raise ValueError(f"{path.name}: expected {columns} columns, found {len(row)}")
            yield [x.strip() for x in row]


def numeric_status(value, unit, known_units):
    if not value:
        return "value_missing"
    try:
        number = float(value)
    except ValueError:
        return "nonnumeric"
    if not math.isfinite(number):
        return "nonfinite"
    if number < 0:
        return "negative"
    if number == 0:
        return "zero"
    if not unit:
        return "positive_unit_missing"
    if unit.lower() not in known_units:
        return "positive_unrecognized_unit"
    return "positive_recognized_unit"


def converted(value, unit, mapping):
    try:
        number = float(value)
    except ValueError:
        return None
    return (
        number * mapping[unit.lower()]
        if math.isfinite(number) and unit.lower() in mapping
        else None
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument(
        "--stage5",
        type=Path,
        required=True,
        help="Directory containing the locked Stage 5 script and aggregate CSVs",
    )
    parser.add_argument("--canada", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    stage5 = args.stage5 or args.root / "03_确认后补充分析" / "阶段五_主要分析"
    if not stage5.exists():
        raise ValueError("Project root is incorrect; pass --root explicitly")
    spec = importlib.util.spec_from_file_location(
        "s5_locked", stage5 / "run_stage5_primary_analysis.py"
    )
    s5 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(s5)
    with (stage5 / "阶段五_主要结局ROR_Fisher_Holm.csv").open(encoding="utf-8-sig", newline="") as f:
        locked = {r["drug"]: r for r in csv.DictReader(f) if r["database"] == "Canada Vigilance"}
    source_files = ["reports.txt", "report_drug.txt", "reactions.txt", "report_links.txt"]
    sources = [
        {
            "file": x,
            "bytes": (args.canada / x).stat().st_size,
            "modified": datetime.fromtimestamp((args.canada / x).stat().st_mtime).isoformat(
                timespec="seconds"
            ),
        }
        for x in source_files
    ]
    with (stage5 / "阶段五_源文件指纹.csv").open(encoding="utf-8-sig", newline="") as f:
        old_sources = {str(Path(r["path"]).resolve()).lower(): r for r in csv.DictReader(f)}
    for row in sources:
        old = old_sources.get(str((args.canada / row["file"]).resolve()).lower())
        row["matches_stage5_size_mtime"] = (
            (row["bytes"] == int(old["bytes"]) and row["modified"] == old["modified"])
            if old
            else None
        )
        if old and not row["matches_stage5_size_mtime"]:
            raise ValueError("Previously used source changed")
    plan = {
        "study_unit": "report, not patient",
        "date": datetime.now().isoformat(timespec="seconds"),
        "source_version": "cvponline_extract_20260228",
        "design": "Descriptive audit of existing spontaneous reports; no follow-up cohort or causal model.",
        "target_drugs": list(s5.DRUG_TERMS),
        "drug_role": "native Suspect",
        "primary_outcome": s5.OUTCOMES["primary_union"],
        "new_tests": 0,
        "raw_records_modified": False,
        "duplicate_links_removed": False,
        "missingness_rule": "Blank values are missing; recorded nonnumeric/nonfinite/zero/negative/unit-missing/unrecognized-unit values are separate, never assumed normal.",
        "duration_rule": "Therapy duration and reaction duration are described separately and never interpreted as time to onset.",
        "body_measure_rule": "The official glossary defines height/weight as recorded when the adverse effect occurred; no pretreatment measurement is established, so no baseline BMI is constructed.",
        "extreme_value_screen": "Raw finite minimum and maximum by field/unit are retained. Body-measure audit flags weight >500 kg or height >300 cm after standard unit conversion as broad plausibility concerns; flags do not remove reports or establish clinical validity for other positive values.",
        "duplicate_link_rule": "Exact REPORT_LINK_NO to REPORT_NO first; exact REPORT_ID mapping is additionally audited, without inferring undocumented normalization. No duplicate report removal or ROR recomputation.",
        "schema_url": SCHEMA_URL,
        "schema_access_date": "2026-10-20",
        "analysis_scope": "Clinical-field availability audit; no new inferential modeling.",
    }
    (args.out / "audit_plan_before_run.json").write_text(
        json.dumps(plan, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("Auditing report keys and clinical-field presence", flush=True)
    reports = {}
    no_to_id = {}
    report_duplicates = Counter()
    body = {}
    raw_missing = defaultdict(Counter)
    units = defaultdict(Counter)
    ranges = {}

    def observe_range(field, value, unit):
        try:
            number = float(value)
        except ValueError:
            return
        if not math.isfinite(number):
            return
        key = (field, unit or "[blank]")
        if key not in ranges:
            ranges[key] = [number, number, 1]
        else:
            r = ranges[key]
            r[0], r[1], r[2] = min(r[0], number), max(r[1], number), r[2] + 1

    file_audit = []
    n = 0
    for row in read_rows(args.canada / "reports.txt", 42):
        n += 1
        rid, no = row[0], row[1]
        report_duplicates["REPORT_ID"] += rid in reports
        report_duplicates["REPORT_NO"] += no in no_to_id
        reports[rid] = no
        no_to_id[no] = rid
        ws = numeric_status(row[19], row[20], WEIGHT_UNITS)
        hs = numeric_status(row[22], row[23], HEIGHT_UNITS)
        raw_missing["weight"][ws] += 1
        raw_missing["height"][hs] += 1
        units["weight"][row[20] or "[blank]"] += 1
        units["height"][row[23] or "[blank]"] += 1
        observe_range("weight", row[19], row[20])
        observe_range("height", row[22], row[23])
        body[rid] = (
            ws,
            hs,
            converted(row[19], row[20], WEIGHT_TO_KG),
            converted(row[22], row[23], HEIGHT_TO_CM),
        )
    if report_duplicates["REPORT_ID"] or report_duplicates["REPORT_NO"]:
        raise ValueError(
            "Duplicate report keys; cannot inherit retained report set without resolution"
        )
    file_audit.append(
        {
            "file": "reports.txt",
            "columns": 42,
            "rows": n,
            "key": "REPORT_ID; REPORT_NO",
            "distinct_primary_keys": len(reports),
            "duplicate_primary_key_rows": report_duplicates["REPORT_ID"],
            "duplicate_report_numbers": report_duplicates["REPORT_NO"],
            "orphan_report_rows": 0,
        }
    )
    target = {d: set() for d in s5.DRUG_TERMS}
    therapy = {d: defaultdict(Counter) for d in target}
    therapy_raw = {d: defaultdict(set) for d in target}
    row_keys = set()
    duplicates = 0
    orphan = 0
    n = 0
    roles = Counter()
    for row in read_rows(args.canada / "report_drug.txt", 22):
        n += 1
        key, rid = row[0], row[1]
        duplicates += key in row_keys
        row_keys.add(key)
        orphan += rid not in reports
        roles[row[4]] += 1
        status = numeric_status(row[17], row[18], DURATION_UNITS)
        raw_missing["therapy_duration"][status] += 1
        units["therapy_duration"][row[18] or "[blank]"] += 1
        observe_range("therapy_duration", row[17], row[18])
        if row[4].upper() != "SUSPECT" or rid not in reports:
            continue
        drug = s5.match_drug(row[3])
        if drug:
            target[drug].add(rid)
            therapy[drug][rid][status] += 1
            if row[17]:
                therapy_raw[drug][rid].add((row[17], row[18]))
    file_audit.append(
        {
            "file": "report_drug.txt",
            "columns": 22,
            "rows": n,
            "key": "REPORT_DRUG_ID",
            "distinct_primary_keys": len(row_keys),
            "duplicate_primary_key_rows": duplicates,
            "orphan_report_rows": orphan,
        }
    )
    del row_keys
    gc.collect()
    print("Auditing primary-reaction durations and locked cohort reconciliation", flush=True)
    primary_terms = {x.lower() for x in s5.OUTCOMES["primary_union"]}
    primary = set()
    primary_duration = defaultdict(Counter)
    primary_duration_raw = defaultdict(set)
    all_target = set.union(*target.values())
    row_keys = set()
    n = duplicates = orphan = 0
    for row in read_rows(args.canada / "reactions.txt", 10):
        n += 1
        key, rid = row[0], row[1]
        duplicates += key in row_keys
        row_keys.add(key)
        orphan += rid not in reports
        status = numeric_status(row[2], row[3], DURATION_UNITS)
        raw_missing["reaction_duration"][status] += 1
        units["reaction_duration"][row[3] or "[blank]"] += 1
        observe_range("reaction_duration", row[2], row[3])
        if rid in reports and row[5].lower() in primary_terms:
            primary.add(rid)
            if rid in all_target:
                primary_duration[rid][status] += 1
                if row[2]:
                    primary_duration_raw[rid].add((row[5], row[2], row[3]))
    file_audit.append(
        {
            "file": "reactions.txt",
            "columns": 10,
            "rows": n,
            "key": "REACTION_ID",
            "distinct_primary_keys": len(row_keys),
            "duplicate_primary_key_rows": duplicates,
            "orphan_report_rows": orphan,
        }
    )
    del row_keys
    gc.collect()
    reconciliation = []
    for d, ids in target.items():
        old = locked[d]
        check = (
            len(reports) == int(old["background_reports_n"])
            and len(ids) == int(old["target_reports_n"])
            and len(primary) == int(old["total_case_reports_n"])
            and len(ids & primary) == int(old["a_target_case_reports"])
        )
        reconciliation.append(
            {
                "drug": d,
                "all_target_reports_n": len(ids),
                "primary_reports_n": len(ids & primary),
                "background_reports_n": len(reports),
                "background_primary_n": len(primary),
                "all_four_locked_counts_match": check,
            }
        )
        if not check:
            raise ValueError(f"Locked Stage 5 count mismatch for {d}")
    save_csv(args.out / "locked_cohort_reconciliation.csv", reconciliation)
    print("Auditing official Duplicate/Linked report associations", flush=True)
    link_types = Counter()
    link_shapes = Counter()
    mapping = Counter()
    typed_links = defaultdict(list)
    row_keys = set()
    semantic_keys = set()
    n = duplicates = orphan = semantic_duplicates = 0
    for row in read_rows(args.canada / "report_links.txt", 5):
        n += 1
        key, rid, typ, _, linked_no = row
        duplicates += key in row_keys
        row_keys.add(key)
        orphan += rid not in reports
        semantic_key = (rid, typ, linked_no)
        semantic_duplicates += semantic_key in semantic_keys
        semantic_keys.add(semantic_key)
        link_types[typ] += 1
        link_shapes[(typ, re.sub("[0-9]", "#", re.sub("[A-Za-z]", "A", linked_no)))] += 1
        by_no = no_to_id.get(linked_no)
        by_id = linked_no if linked_no in reports else None
        if by_no and by_id and by_no != by_id:
            mode = "ambiguous_REPORT_NO_vs_REPORT_ID"
            endpoint = None
        elif by_no:
            mode = "exact_REPORT_NO"
            endpoint = by_no
        elif by_id:
            mode = "exact_REPORT_ID_only"
            endpoint = by_id
        else:
            mode = "no_exact_retained_report_match" if linked_no else "empty_link_number"
            endpoint = None
        mapping[(typ, mode)] += 1
        typed_links[typ].append((rid, endpoint, mode, linked_no))
    file_audit.append(
        {
            "file": "report_links.txt",
            "columns": 5,
            "rows": n,
            "key": "REPORT_LINK_ID",
            "distinct_primary_keys": len(row_keys),
            "duplicate_primary_key_rows": duplicates,
            "duplicate_report_type_link_number_rows": semantic_duplicates,
            "orphan_report_rows": orphan,
        }
    )
    del row_keys, semantic_keys
    gc.collect()
    # The file/schema/cohort audit is recorded before any aggregate interpretation.
    audit = {
        "status": "schema_keys_and_locked_cohort_checked",
        "files": file_audit,
        "sources": sources,
        "drug_role_counts": dict(roles),
        "link_type_counts": dict(link_types),
        "all_locked_cohort_counts_match": all(
            r["all_four_locked_counts_match"] for r in reconciliation
        ),
        "no_patient_records_exported": True,
    }
    (args.out / "source_schema_key_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    save_csv(
        args.out / "raw_field_status_counts.csv",
        [
            {
                "field": field,
                "status": status,
                "row_n": count,
                "field_table_rows_n": sum(c.values()),
                "percent": round(100 * count / sum(c.values()), 4),
            }
            for field, c in raw_missing.items()
            for status, count in c.items()
        ],
    )
    save_csv(
        args.out / "raw_unit_counts.csv",
        [
            {"field": field, "unit": unit, "row_n": count}
            for field, c in units.items()
            for unit, count in c.items()
        ],
    )
    save_csv(
        args.out / "raw_numeric_ranges_by_unit.csv",
        [
            {
                "field": field,
                "unit": unit,
                "finite_value_rows_n": v[2],
                "minimum": v[0],
                "maximum": v[1],
                "no_unit_conversion_for_duration": True,
            }
            for (field, unit), v in ranges.items()
        ],
    )
    save_csv(
        args.out / "link_mapping_counts.csv",
        [
            {"record_type": typ, "mapping": mode, "rows_n": count}
            for (typ, mode), count in mapping.items()
        ],
    )
    save_csv(
        args.out / "link_number_format_counts.csv",
        [
            {"record_type": typ, "format_only_no_identifiers": shape, "rows_n": count}
            for (typ, shape), count in link_shapes.items()
        ],
    )
    summary = []
    patterns = []
    link_summary = []
    graph_summary = []
    status_order = [
        "value_missing",
        "nonnumeric",
        "nonfinite",
        "negative",
        "zero",
        "positive_unit_missing",
        "positive_unrecognized_unit",
        "positive_recognized_unit",
    ]
    for d, all_ids in target.items():
        for subset, ids in [
            ("all_target_reports", all_ids),
            ("primary_outcome_reports", all_ids & primary),
        ]:
            den = len(ids)
            ws = Counter(body[rid][0] for rid in ids)
            hs = Counter(body[rid][1] for rid in ids)
            for field, counts in [("weight", ws), ("height", hs)]:
                for status in status_order:
                    count = counts[status]
                    summary.append(
                        {
                            "drug": d,
                            "subset": subset,
                            "variable": field,
                            "measure": status,
                            "reports_n": count,
                            "denominator_reports_n": den,
                            "percent": round(100 * count / den, 4) if den else "",
                        }
                    )
            both = sum(
                body[rid][:2] == ("positive_recognized_unit", "positive_recognized_unit")
                for rid in ids
            )
            summary.append(
                {
                    "drug": d,
                    "subset": subset,
                    "variable": "height_and_weight",
                    "measure": "both_positive_recognized_unit",
                    "reports_n": both,
                    "denominator_reports_n": den,
                    "percent": round(100 * both / den, 4) if den else "",
                }
            )
            for field, index, limit in [("weight", 2, 500), ("height", 3, 300)]:
                flag = sum(body[rid][index] is not None and body[rid][index] > limit for rid in ids)
                summary.append(
                    {
                        "drug": d,
                        "subset": subset,
                        "variable": field,
                        "measure": "above_broad_plausibility_screen",
                        "reports_n": flag,
                        "denominator_reports_n": den,
                        "percent": round(100 * flag / den, 4) if den else "",
                    }
                )
            therapy_pattern = Counter()
            for rid in ids:
                c = therapy[d][rid]
                raw_present = any(k != "value_missing" and v for k, v in c.items())
                usable = bool(c["positive_recognized_unit"])
                therapy_pattern[
                    (
                        body[rid][0] == "value_missing",
                        body[rid][1] == "value_missing",
                        not raw_present,
                        not usable,
                    )
                ] += 1
            for (wm, hm, tm, tu), count in therapy_pattern.items():
                patterns.append(
                    {
                        "drug": d,
                        "subset": subset,
                        "weight_value_missing": wm,
                        "height_value_missing": hm,
                        "all_target_suspect_therapy_values_missing": tm,
                        "no_positive_known_unit_target_suspect_therapy": tu,
                        "reports_n": count,
                        "denominator_reports_n": den,
                    }
                )
            for field, records, raws in [
                ("target_suspect_therapy_duration", therapy[d], therapy_raw[d]),
                ("primary_reaction_duration", primary_duration, primary_duration_raw),
            ]:
                eligible = ids if field.startswith("target") else ids & primary
                ed = len(eligible)
                has_any = sum(
                    any(k != "value_missing" and v for k, v in records[rid].items())
                    for rid in eligible
                )
                has_good = sum(bool(records[rid]["positive_recognized_unit"]) for rid in eligible)
                multi = sum(len(raws[rid]) > 1 for rid in eligible)
                measures = {
                    "any_value_recorded": has_any,
                    "all_values_missing": ed - has_any,
                    "any_positive_recognized_unit": has_good,
                    "no_positive_recognized_unit": ed - has_good,
                    "multiple_distinct_recorded_value_unit_combinations": multi,
                }
                for key, count in measures.items():
                    summary.append(
                        {
                            "drug": d,
                            "subset": subset,
                            "variable": field,
                            "measure": key,
                            "reports_n": count,
                            "denominator_reports_n": ed,
                            "percent": round(100 * count / ed, 4) if ed else "",
                        }
                    )
            for typ, links in typed_links.items():
                outgoing = [(a, b, m, no) for a, b, m, no in links if a in ids]
                incoming = [(a, b, m, no) for a, b, m, no in links if b in ids]
                involved = {a for a, b, m, no in outgoing} | {b for a, b, m, no in incoming}
                relevant_edges = outgoing + incoming
                pairs = {
                    tuple(sorted([a, b]))
                    for a, b, m, no in relevant_edges
                    if a in reports and b in reports and a != b
                }
                outgoing_unmatched = {a for a, b, m, no in outgoing if b is None}
                both_cohort_pairs = {pair for pair in pairs if all(x in ids for x in pair)}
                same_drug_pairs = {pair for pair in pairs if all(x in all_ids for x in pair)}
                primary_same_drug_pairs = {
                    pair for pair in pairs if all(x in (all_ids & primary) for x in pair)
                }
                link_summary.append(
                    {
                        "drug": d,
                        "subset": subset,
                        "record_type": typ,
                        "denominator_reports_n": den,
                        "outgoing_link_rows_n": len(outgoing),
                        "outgoing_linked_source_reports_n": len({a for a, b, m, no in outgoing}),
                        "incoming_link_rows_n": len(incoming),
                        "reports_involved_either_direction_n": len(involved),
                        "reports_involved_percent": round(100 * len(involved) / den, 4)
                        if den
                        else "",
                        "outgoing_source_reports_with_unmatched_endpoint_n": len(
                            outgoing_unmatched
                        ),
                        "outgoing_rows_unmatched_endpoint_n": sum(
                            b is None for a, b, m, no in outgoing
                        ),
                        "retained_distinct_undirected_pairs_touching_subset_n": len(pairs),
                        "retained_distinct_pairs_both_in_subset_n": len(both_cohort_pairs),
                        "retained_distinct_pairs_both_same_target_drug_n": len(same_drug_pairs),
                        "retained_distinct_pairs_both_same_target_drug_primary_n": len(
                            primary_same_drug_pairs
                        ),
                        "self_link_rows_n": sum(a == b for a, b, m, no in outgoing),
                        "ambiguous_link_number_mapping_rows_n": sum(
                            m == "ambiguous_REPORT_NO_vs_REPORT_ID" for a, b, m, no in outgoing
                        ),
                    }
                )
    # Global Duplicate connected components are descriptive; they are not adjudicated patients.
    duplicate_links = [
        x for typ, links in typed_links.items() if typ.lower() == "duplicate" for x in links
    ]
    edges = {
        tuple(sorted((a, b)))
        for a, b, m, no in duplicate_links
        if a in reports and b in reports and a != b
    }
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        aa, bb = find(a), find(b)
        if aa != bb:
            parent[bb] = aa
    components = defaultdict(set)
    for rid in parent:
        components[find(rid)].add(rid)
    for d, ids in target.items():
        for subset, cohort in [
            ("all_target_reports", ids),
            ("primary_outcome_reports", ids & primary),
        ]:
            touching = [group for group in components.values() if group & cohort]
            multi_same = [group for group in touching if len(group & cohort) > 1]
            graph_summary.append(
                {
                    "drug": d,
                    "subset": subset,
                    "denominator_reports_n": len(cohort),
                    "duplicate_components_touching_subset_n": len(touching),
                    "reports_in_duplicate_components_n": sum(
                        len(group & cohort) for group in touching
                    ),
                    "components_with_multiple_subset_reports_n": len(multi_same),
                    "possible_redundant_subset_reports_if_one_per_component": sum(
                        len(group & cohort) - 1 for group in multi_same
                    ),
                    "no_report_removed": True,
                }
            )
    save_csv(args.out / "clinical_field_availability_by_drug.csv", summary)
    save_csv(args.out / "clinical_field_missingness_patterns.csv", patterns)
    save_csv(args.out / "duplicate_and_linked_report_summary.csv", link_summary)
    save_csv(args.out / "duplicate_component_summary.csv", graph_summary)
    validation = {
        "status": "completed_aggregate_audit",
        "all_locked_counts_match": True,
        "report_key_duplicates": dict(report_duplicates),
        "new_models_or_hypothesis_tests": 0,
        "report_removals": 0,
        "duplicate_link_rows": len(duplicate_links),
        "retained_duplicate_distinct_undirected_edges": len(edges),
        "retained_duplicate_components": len(components),
        "duplicate_component_size_distribution": dict(Counter(len(v) for v in components.values())),
        "limitations": [
            "A positive numeric value with a recognized unit is a syntactic availability check, not clinical validation.",
            "The official glossary defines height/weight as recorded when the adverse effect occurred, not as pretreatment measurements; baseline BMI and nutritional status cannot be validated.",
            "Therapy duration and reaction duration are different fields and neither is an onset latency.",
            "Duplicate links identify source-labelled associations, not verified patients; unresolved endpoints do not establish absence from the source agency system.",
            "No link-based report removal, disproportionality recomputation, or imputation was performed.",
            "Co-recorded information cannot establish pretreatment status or temporal causal ordering.",
        ],
        "dependency_script_sha256": hashlib.sha256(
            (stage5 / "run_stage5_primary_analysis.py").read_bytes()
        ).hexdigest(),
    }
    (args.out / "audit_validation.json").write_text(
        json.dumps(validation, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "status": validation["status"],
                "source_rows": {x["file"]: x["rows"] for x in file_audit},
                "cohort": reconciliation,
                "duplicate_summary": graph_summary,
            },
            ensure_ascii=False,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
