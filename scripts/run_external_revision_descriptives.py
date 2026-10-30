"""Audited report-level external-database descriptions and four-domain overlap.

Reads source files only. Writes aggregates; no report identifiers are exported.
Run with Python 3.10+; inherits locked Stage 5 drug/PT definitions and counts.
"""
from __future__ import annotations
import argparse, csv, gc, hashlib, importlib.util, itertools, json, math, re
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from package_paths import OUTPUT, CANADA, JADER

SCHEMA_URL = "https://www.canada.ca/en/health-canada/services/drugs-health-products/medeffect-canada/adverse-reaction-database/canada-vigilance-adverse-reaction-online-database-data-structures.html"
DOMAINS = {
    "primary": ["Diabetic ketoacidosis", "Ketoacidosis", "Metabolic acidosis"],
    "dehydration": ["Dehydration", "Hypovolaemia"],
    "electrolyte": ["Hypokalaemia", "Hyponatraemia"],
    "renal": ["Acute kidney injury", "Renal failure", "Renal impairment"],
}
CAN_AGE_GROUPS = [
    "<18 years",
    "18-39 years",
    "40-64 years",
    "65-74 years",
    ">=75 years",
    "Missing",
    "Invalid",
]
UNKNOWN_INDICATIONS = {
    "unknown",
    "unknown indication",
    "not reported",
    "product used for unknown indication",
    "drug use for unknown indication",
    "不明",
    "適応症不明",
    "使用理由不明",
}
JP_AGES = {
    "10歳未満": "<10 years",
    **{f"{n}歳代": f"{n}-{n+9} years" for n in range(10, 110, 10)},
    "高齢者": "Elderly (no decade)",
    "成人": "Adult (no decade)",
    "小児": "Child (no decade)",
    "青少年": "Adolescent (no decade)",
    "新生児": "Neonate",
    "乳幼児": "Infant/toddler",
    "胎児": "Fetus",
    "第一トライメスター": "First trimester",
    "第二トライメスター": "Second trimester",
    "第三トライメスター": "Third trimester",
}


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def save_csv(path, rows):
    if not rows:
        return
    columns = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def pct(n, d):
    return round(100 * n / d, 4) if d else ""


def finite_number(value):
    try:
        n = float(value)
        return n if math.isfinite(n) else None
    except (TypeError, ValueError):
        return None


def age_category(db, row):
    if db == "Canada Vigilance":
        raw = row[13].strip()  # official AGE_Y; do not treat AGE as years
        if not raw:
            return "Missing", "missing", ""
        value = finite_number(raw)
        if value is None or not 0 <= value <= 120:
            return "Invalid", "invalid", raw
        group = (
            "<18 years"
            if value < 18
            else "18-39 years"
            if value < 40
            else "40-64 years"
            if value < 65
            else "65-74 years"
            if value < 75
            else ">=75 years"
        )
        return group, "numeric_years", raw
    raw = row[3].strip()
    if raw in {"", "不明"}:
        return "Missing" if not raw else "Unknown", "missing" if not raw else "unknown", raw
    if raw in JP_AGES:
        return (
            JP_AGES[raw],
            "age_interval" if ("歳代" in raw or raw == "10歳未満") else "noninterval_category",
            raw,
        )
    return "Unrecognized age category", "invalid", raw


def sex_category(db, row, gender_codes):
    raw = row[10].strip() if db == "Canada Vigilance" else row[2].strip()
    if not raw:
        return "Missing", "missing", raw
    female, male = ("Female", "Male") if db == "Canada Vigilance" else ("女性", "男性")
    if raw == female:
        return "Female", "recorded", raw
    if raw == male:
        return "Male", "recorded", raw
    if raw in {"Unknown", "不明"}:
        return "Unknown", "unknown", raw
    if db == "Canada Vigilance" and raw in gender_codes.values():
        return raw, "other_recorded", raw
    return "Unrecognized sex code", "invalid", raw


def load_canada(s5, root):
    reports = {}
    report_numbers = set()
    report_lengths = Counter()
    duplicate_ids = 0
    duplicate_numbers = 0
    with (root / "reports.txt").open(encoding="utf-8-sig", newline="") as f:
        for row in csv.reader(f, delimiter="$"):
            if not row:
                continue
            report_lengths[len(row)] += 1
            if len(row) != 42:
                raise ValueError("Canada reports schema not 42 columns")
            rid = row[0].strip()
            duplicate_ids += rid in reports
            duplicate_numbers += row[1] in report_numbers
            if rid in reports:
                raise ValueError("Canada duplicate report_id; stop before choosing a version")
            reports[rid] = (
                row[2],
            )  # demographics are reread only for the final small primary cohort
            report_numbers.add(row[1])
    del report_numbers
    target = {d: set() for d in s5.DRUG_TERMS}
    links = {}
    roles = Counter()
    drug_rows = 0
    blank_rows = Counter()
    with (root / "report_drug.txt").open(encoding="utf-8-sig", newline="") as f:
        for row in csv.reader(f, delimiter="$"):
            if not row:
                blank_rows["report_drug.txt"] += 1
                continue
            drug_rows += 1
            if len(row) != 22:
                raise ValueError("Canada report_drug schema not 22 columns")
            roles[row[4].strip()] += 1
            if row[4].strip().upper() != "SUSPECT" or row[1].strip() not in reports:
                continue
            drug = s5.match_drug(row[3])
            if drug:
                rid = row[1].strip()
                target[drug].add(rid)
                key = row[0].strip()
                value = (rid, row[2].strip(), drug)
                if key in links and links[key] != value:
                    raise ValueError("Ambiguous Canada REPORT_DRUG_ID")
                links[key] = value
    pts = {x.lower(): set() for terms in DOMAINS.values() for x in terms}
    reaction_rows = 0
    orphan_reactions = 0
    with (root / "reactions.txt").open(encoding="utf-8-sig", newline="") as f:
        for row in csv.reader(f, delimiter="$"):
            if not row:
                blank_rows["reactions.txt"] += 1
                continue
            reaction_rows += 1
            if len(row) != 10:
                raise ValueError("Canada reactions schema not 10 columns")
            rid, pt = row[1].strip(), row[5].strip().lower()
            if rid not in reports:
                orphan_reactions += 1
                continue
            if pt in pts:
                pts[pt].add(rid)
    primary = set.union(*(pts[p.lower()] for p in DOMAINS["primary"]))
    primary_targets = set.union(*(ids & primary for ids in target.values()))
    demographic = {}
    units = Counter()
    gender_mismatches = 0
    age_year_unit_mismatches = 0
    gender_codes = {}
    with (root / "gender_lx.txt").open(encoding="utf-8-sig", newline="") as f:
        for row in csv.reader(f, delimiter="$"):
            if row:
                gender_codes[row[1].strip()] = row[2].strip()
    with (root / "reports.txt").open(encoding="utf-8-sig", newline="") as f:
        for row in csv.reader(f, delimiter="$"):
            if not row or row[0].strip() not in primary_targets:
                continue
            rid = row[0].strip()
            demographic[rid] = row
            units[row[14].strip() or "[blank]"] += 1
            if row[9].strip() and gender_codes.get(row[9].strip()) != row[10].strip():
                gender_mismatches += 1
            age, age_y = finite_number(row[12]), finite_number(row[13])
            if (
                row[14].strip() == "Years"
                and age is not None
                and age_y is not None
                and abs(age - age_y) > 0.001
            ):
                age_year_unit_mismatches += 1
    indications = {d: defaultdict(set) for d in target}
    indication_rows = 0
    matched_rows = 0
    indication_mismatches = 0
    with (root / "report_drug_indication.txt").open(encoding="utf-8-sig", newline="") as f:
        for row in csv.reader(f, delimiter="$"):
            if not row:
                blank_rows["report_drug_indication.txt"] += 1
                continue
            indication_rows += 1
            if len(row) != 6:
                raise ValueError("Canada indication schema not 6 columns")
            matched = links.get(row[0].strip())
            if not matched:
                continue
            rid, product_id, drug = matched
            if row[1].strip() != rid or row[2].strip() != product_id:
                indication_mismatches += 1
                continue
            if rid in primary:
                matched_rows += 1
                if row[4].strip():
                    indications[drug][rid].add(row[4].strip())
    if indication_mismatches or gender_mismatches or age_year_unit_mismatches:
        raise ValueError("Canada target data linkage or field validation failed")
    meta = {
        "database": "Canada Vigilance",
        "source_version": "cvponline_extract_20260228",
        "reports_n": len(reports),
        "duplicate_report_ids": duplicate_ids,
        "duplicate_report_numbers": duplicate_numbers,
        "report_row_lengths": dict(report_lengths),
        "drug_rows": drug_rows,
        "blank_lines_skipped": dict(blank_rows),
        "drug_role_counts": dict(roles),
        "reaction_rows": reaction_rows,
        "orphan_reaction_rows": orphan_reactions,
        "primary_target_age_units": dict(units),
        "gender_lookup": gender_codes,
        "gender_code_text_mismatches": gender_mismatches,
        "years_age_vs_age_y_mismatches": age_year_unit_mismatches,
        "indication_rows": indication_rows,
        "target_primary_matched_indication_rows": matched_rows,
        "indication_exact_link_mismatches": indication_mismatches,
        "indication_link": "REPORT_DRUG_ID plus REPORT_ID and DRUG_PRODUCT_ID, native Suspect role",
        "deduplication": "One row per REPORT_ID and REPORT_NO verified; report-level sets. Cross-ID clinical duplicates not resolved.",
    }
    return len(reports), target, pts, demographic, indications, meta, gender_codes


def load_jader(s5, root):
    latest = {}
    counts = Counter()
    demo_rows = 0
    duplicate_keys = 0
    seen = set()
    with (root / "demo202605.csv").open(encoding="cp932", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        if header[:4] != ["識別番号", "報告回数", "性別", "年齢"]:
            raise ValueError("JADER demographic schema changed")
        for row in reader:
            demo_rows += 1
            counts[len(row)] += 1
            key = s5.jkey(row[0], row[1])
            duplicate_keys += key in seen
            seen.add(key)
            version = int(row[1])
            if row[0] not in latest or version > latest[row[0]]:
                latest[row[0]] = version
    if duplicate_keys:
        raise ValueError("Duplicate JADER ID/version demographic rows")
    reports = {s5.jkey(rid, version) for rid, version in latest.items()}
    del latest, seen
    target = {d: set() for d in s5.DRUG_TERMS}
    indication_all = {d: defaultdict(set) for d in target}
    roles = Counter()
    drug_rows = 0
    with (root / "drug202605.csv").open(encoding="cp932", newline="") as f:
        reader = csv.reader(f)
        drug_header = next(reader)
        if drug_header[12] != "使用理由":
            raise ValueError("JADER indication column changed")
        for row in reader:
            drug_rows += 1
            roles[row[3].strip()] += 1
            key = s5.jkey(row[0], row[1])
            if row[3].strip() != "被疑薬" or key not in reports:
                continue
            drug = s5.match_drug(row[4], row[5])
            if drug:
                target[drug].add(key)
                if row[12].strip():
                    indication_all[drug][key].add(row[12].strip())
    pts = {x.lower(): set() for terms in DOMAINS.values() for x in terms}
    reaction_rows = 0
    orphan_reactions = 0
    with (root / "reac202605.csv").open(encoding="cp932", newline="") as f:
        reader = csv.reader(f)
        reaction_header = next(reader)
        for row in reader:
            reaction_rows += 1
            key = s5.jkey(row[0], row[1])
            if key not in reports:
                orphan_reactions += 1
                continue
            pt = s5.JP_TO_EN.get(row[3].strip(), "").lower()
            if pt in pts:
                pts[pt].add(key)
    primary = set.union(*(pts[p.lower()] for p in DOMAINS["primary"]))
    primary_targets = set.union(*(ids & primary for ids in target.values()))
    demographic = {}
    with (root / "demo202605.csv").open(encoding="cp932", newline="") as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            key = s5.jkey(row[0], row[1])
            if key in primary_targets:
                demographic[key] = row
    indications = {
        d: {rid: terms for rid, terms in rows.items() if rid in primary}
        for d, rows in indication_all.items()
    }
    meta = {
        "database": "JADER",
        "source_version": "202605",
        "reports_n": len(reports),
        "demo_rows": demo_rows,
        "superseded_demo_rows": demo_rows - len(reports),
        "duplicate_id_version_rows": duplicate_keys,
        "demo_row_lengths": dict(counts),
        "demo_header": header,
        "drug_header": drug_header,
        "reaction_header": reaction_header,
        "drug_rows": drug_rows,
        "drug_role_counts": dict(roles),
        "reaction_rows": reaction_rows,
        "orphan_reaction_rows": orphan_reactions,
        "indication_link": "使用理由 on the exact target-suspect drug row, matched REPORT ID and latest version",
        "deduplication": "Highest numeric report version for each ID; report-level sets. No clinical cross-ID duplicate inference.",
    }
    return len(reports), target, pts, demographic, indications, meta, {}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=OUTPUT)
    parser.add_argument("--canada", type=Path, default=CANADA)
    parser.add_argument("--jader", type=Path, default=JADER)
    parser.add_argument(
        "--out", type=Path, default=OUTPUT / "03_确认后补充分析" / "阶段八_核验后修复" / "External"
    )
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    stage5 = args.root / "03_确认后补充分析" / "阶段五_主要分析"
    spec = importlib.util.spec_from_file_location(
        "locked_stage5_external", Path(__file__).resolve().parent / "run_stage5_primary_analysis.py"
    )
    s5 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(s5)
    locked = {
        (r["database"], r["drug"]): r for r in read_csv(stage5 / "阶段五_主要结局ROR_Fisher_Holm.csv")
    }
    secondary = {
        (r["database"], r["drug"], r["outcome"]): r
        for r in read_csv(stage5 / "阶段五_次要结局ROR_Fisher_BY.csv")
    }
    source_locked = {
        str(Path(r["path"]).resolve()).lower(): r for r in read_csv(stage5 / "阶段五_源文件指纹.csv")
    }
    sources = []
    audit = []
    reconcile = []
    demographic_summary = []
    distribution = []
    indication_rows = []
    domain_rows = []
    combination_rows = []
    intersection_rows = []
    broad_rows = []
    missing_patterns = []
    for db, folder, loader, files in [
        (
            "Canada Vigilance",
            args.canada,
            load_canada,
            [
                "reports.txt",
                "report_drug.txt",
                "reactions.txt",
                "report_drug_indication.txt",
                "gender_lx.txt",
            ],
        ),
        ("JADER", args.jader, load_jader, ["demo202605.csv", "drug202605.csv", "reac202605.csv"]),
    ]:
        print(f"Auditing {db} sources and report counts", flush=True)
        for filename in files:
            p = folder / filename
            stat = p.stat()
            modified = datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds")
            old = source_locked.get(str(p.resolve()).lower())
            same = (
                stat.st_size == int(old["bytes"]) and modified == old["modified"] if old else None
            )
            if old and not same:
                raise ValueError(f"Source fingerprint changed: {filename}")
            sources.append(
                {
                    "database": db,
                    "file": str(p),
                    "bytes": stat.st_size,
                    "modified": modified,
                    "matches_stage5_size_mtime": same,
                    "fingerprint_type": "size+mtime; not cryptographic raw-data checksum",
                }
            )
        N, target, pts, demographic, indications, meta, gender_codes = loader(s5, folder)
        domain_sets = {
            key: set.union(*(pts[p.lower()] for p in terms)) for key, terms in DOMAINS.items()
        }
        broad = set.union(*domain_sets.values())
        for drug, exposed in target.items():
            prior = locked[(db, drug)]
            primary = exposed & domain_sets["primary"]
            old_broad = secondary[(db, drug, "severe_metabolic_renal_composite")]
            checks = [
                N == int(prior["background_reports_n"]),
                len(exposed) == int(prior["target_reports_n"]),
                len(primary) == int(prior["a_target_case_reports"]),
                len(domain_sets["primary"]) == int(prior["total_case_reports_n"]),
                len(exposed & broad) == int(old_broad["a_target_case_reports"]),
                len(broad) == int(old_broad["total_case_reports_n"]),
            ]
            if not all(checks):
                raise ValueError(f"Locked Stage 5 mismatch {db}/{drug}")
            reconcile.append(
                {
                    "database": db,
                    "drug": drug,
                    "background_reports_n": N,
                    "target_reports_n": len(exposed),
                    "primary_reports_n": len(primary),
                    "broad_reports_n": len(exposed & broad),
                    "all_six_stage5_checks_pass": all(checks),
                }
            )
        # Only after source/schema/role/duplicate/count audits pass, produce new descriptions.
        meta["stage5_all_drugs_reconcile"] = True
        audit.append(meta)
        save_csv(args.out / "external_stage5_reconciliation.csv", reconcile)
        (args.out / "external_data_audit.json").write_text(
            json.dumps(
                {
                    "status": "audit_passed_before_description",
                    "sources": sources,
                    "databases": audit,
                },
                ensure_ascii=False,
                indent=2,
            ),
            encoding="utf-8",
        )
        print(f"{db} audit passed; computing descriptive summaries", flush=True)
        for drug, exposed in target.items():
            primary = exposed & domain_sets["primary"]
            n = len(primary)
            broad_target = exposed & broad
            bn = len(broad_target)
            age_counts = Counter()
            sex_counts = Counter()
            age_status = Counter()
            sex_status = Counter()
            patterns = Counter()
            ages = []
            terms_counts = Counter()
            indication_known = set()
            indication_unknown = set()
            for rid in primary:
                row = demographic[rid]
                age, astat, araw = age_category(db, row)
                sex, sstat, sraw = sex_category(db, row, gender_codes)
                age_counts[(age, astat, araw if db == "JADER" else "")] += 1
                sex_counts[(sex, sstat, sraw)] += 1
                age_status[astat] += 1
                sex_status[sstat] += 1
                if astat == "numeric_years":
                    ages.append(float(araw))
                terms = indications[drug].get(rid, set())
                known = {t for t in terms if t.lower() not in UNKNOWN_INDICATIONS}
                if known:
                    indication_known.add(rid)
                elif terms:
                    indication_unknown.add(rid)
                for term in terms:
                    terms_counts[term] += 1
                patterns[
                    (
                        astat in {"missing", "unknown", "invalid"},
                        sstat in {"missing", "unknown", "invalid"},
                        not bool(known),
                    )
                ] += 1
            if sum(age_counts.values()) != n or sum(sex_counts.values()) != n:
                raise ValueError("Demographic count reconciliation failed")
            summary = {
                "database": db,
                "drug": drug,
                "primary_reports_n": n,
                "age_valid_numeric_or_interval_n": age_status["numeric_years"]
                + age_status["age_interval"],
                "age_broad_category_only_n": age_status["noninterval_category"],
                "age_missing_or_unknown_n": age_status["missing"] + age_status["unknown"],
                "age_invalid_n": age_status["invalid"],
                "female_n": sum(v for (k, _, _), v in sex_counts.items() if k == "Female"),
                "male_n": sum(v for (k, _, _), v in sex_counts.items() if k == "Male"),
                "sex_unknown_n": sex_status["unknown"],
                "sex_blank_n": sex_status["missing"],
                "sex_other_or_invalid_n": sex_status["other_recorded"] + sex_status["invalid"],
                "indication_recorded_known_n": len(indication_known),
                "indication_recorded_unknown_only_n": len(indication_unknown),
                "indication_absent_n": n - len(indication_known) - len(indication_unknown),
                "indication_missing_or_unknown_pct": pct(n - len(indication_known), n),
            }
            for field in [
                "age_valid_numeric_or_interval_n",
                "age_broad_category_only_n",
                "age_missing_or_unknown_n",
                "age_invalid_n",
                "female_n",
                "male_n",
                "sex_unknown_n",
                "sex_blank_n",
                "sex_other_or_invalid_n",
            ]:
                summary[field[:-2] + "_pct"] = pct(summary[field], n)
            if ages:
                import statistics

                summary.update(
                    age_median_years=statistics.median(ages),
                    age_min_years=min(ages),
                    age_max_years=max(ages),
                )
            else:
                summary.update(age_median_years="", age_min_years="", age_max_years="")
            demographic_summary.append(summary)
            if db == "Canada Vigilance":
                for category in CAN_AGE_GROUPS:
                    status = (
                        "missing"
                        if category == "Missing"
                        else "invalid"
                        if category == "Invalid"
                        else "numeric_years"
                    )
                    age_counts.setdefault((category, status, ""), 0)
            else:
                for raw, category in JP_AGES.items():
                    age_counts.setdefault(
                        (
                            category,
                            "age_interval"
                            if ("歳代" in raw or raw == "10歳未満")
                            else "noninterval_category",
                            raw,
                        ),
                        0,
                    )
                age_counts.setdefault(("Unknown", "unknown", "不明"), 0)
                age_counts.setdefault(("Missing", "missing", ""), 0)
            for category, status in [
                ("Female", "recorded"),
                ("Male", "recorded"),
                ("Unknown", "unknown"),
                ("Missing", "missing"),
            ]:
                raw = (
                    category
                    if db == "Canada Vigilance"
                    else {"Female": "女性", "Male": "男性", "Unknown": "不明", "Missing": ""}[category]
                )
                if category == "Missing":
                    raw = ""
                sex_counts.setdefault((category, status, raw), 0)
            for variable, counts in [("age", age_counts), ("sex", sex_counts)]:
                for (category, status, raw), count in counts.items():
                    distribution.append(
                        {
                            "database": db,
                            "drug": drug,
                            "variable": variable,
                            "category": category,
                            "source_category": raw,
                            "status": status,
                            "reports_n": count,
                            "primary_reports_n": n,
                            "percent_all_primary": pct(count, n),
                        }
                    )
            for term, count in terms_counts.items():
                indication_rows.append(
                    {
                        "database": db,
                        "drug": drug,
                        "target_drug_indication_raw": term,
                        "indication_status": "explicit_unknown"
                        if term.lower() in UNKNOWN_INDICATIONS
                        else "recorded_term_not_clinically_validated",
                        "reports_n": count,
                        "primary_reports_n": n,
                        "percent_all_primary": pct(count, n),
                        "mutually_exclusive": False,
                    }
                )
            indication_rows.append(
                {
                    "database": db,
                    "drug": drug,
                    "target_drug_indication_raw": "[no known indication]",
                    "indication_status": "missing_or_unknown_summary",
                    "reports_n": n - len(indication_known),
                    "primary_reports_n": n,
                    "percent_all_primary": pct(n - len(indication_known), n),
                    "mutually_exclusive": False,
                }
            )
            for (am, sm, im), count in patterns.items():
                missing_patterns.append(
                    {
                        "database": db,
                        "drug": drug,
                        "age_missing_unknown_or_invalid": am,
                        "sex_missing_unknown_or_invalid": sm,
                        "indication_missing_unknown": im,
                        "reports_n": count,
                        "primary_reports_n": n,
                        "percent_all_primary": pct(count, n),
                    }
                )
            within = {domain: ids & exposed for domain, ids in domain_sets.items()}
            names = list(DOMAINS)
            combo_counts = Counter(
                sum(1 << i for i, domain in enumerate(names) if rid in within[domain])
                for rid in broad_target
            )
            if sum(combo_counts.values()) != bn:
                raise ValueError("Mutually exclusive combinations do not sum to broad union")
            for domain, ids in within.items():
                domain_rows.append(
                    {
                        "database": db,
                        "drug": drug,
                        "domain": domain,
                        "reports_n": len(ids),
                        "broad_reports_n": bn,
                        "percent_broad": pct(len(ids), bn),
                        "target_reports_n": len(exposed),
                        "domains_mutually_exclusive": False,
                    }
                )
            for mask in range(1, 16):
                members = [name for i, name in enumerate(names) if mask & (1 << i)]
                combination_rows.append(
                    {
                        "database": db,
                        "drug": drug,
                        "combination_mask": mask,
                        "exact_domains": " + ".join(members),
                        "excluded_domains": " + ".join(
                            name for name in names if name not in members
                        ),
                        "reports_n": combo_counts[mask],
                        "broad_reports_n": bn,
                        "percent_broad": pct(combo_counts[mask], bn),
                        "mutually_exclusive": True,
                    }
                )
            for size in [2, 3, 4]:
                for members in itertools.combinations(names, size):
                    ids = set.intersection(*(within[member] for member in members))
                    mask = sum(1 << names.index(member) for member in members)
                    recovered = sum(
                        count for pattern, count in combo_counts.items() if pattern & mask == mask
                    )
                    if len(ids) != recovered:
                        raise ValueError("Intersection does not match exclusive combinations")
                    intersection_rows.append(
                        {
                            "database": db,
                            "drug": drug,
                            "intersection_order": size,
                            "domains": " & ".join(members),
                            "reports_n": len(ids),
                            "broad_reports_n": bn,
                            "percent_broad": pct(len(ids), bn),
                            "other_domains_allowed": True,
                        }
                    )
            broad_rows.append(
                {
                    "database": db,
                    "drug": drug,
                    "target_reports_n": len(exposed),
                    "broad_reports_n": bn,
                    **{f"{domain}_n": len(ids) for domain, ids in within.items()},
                    "exclusive_combinations_sum_n": sum(combo_counts.values()),
                    "all_intersections_reconstruct": True,
                }
            )
        del target, pts, demographic, indications, domain_sets, broad
        gc.collect()
    outputs = {
        "external_primary_demographic_summary.csv": demographic_summary,
        "external_primary_demographic_distribution.csv": distribution,
        "external_primary_indication_distribution.csv": indication_rows,
        "external_primary_missingness_patterns.csv": missing_patterns,
        "external_composite_domain_contributions.csv": domain_rows,
        "external_composite_15_exclusive_combinations.csv": combination_rows,
        "external_composite_all_intersections.csv": intersection_rows,
        "external_composite_reconciliation.csv": broad_rows,
        "external_source_provenance.csv": sources,
    }
    for filename, rows in outputs.items():
        save_csv(args.out / filename, rows)
    final = {
        "status": "completed_after_audit",
        "date": datetime.now().isoformat(timespec="seconds"),
        "study_unit": "native-suspect report, not patient",
        "purpose": "Revision-stage descriptive analysis; not prespecified at original study inception",
        "new_hypothesis_tests": 0,
        "stage5_drug_database_reconciliations": len(reconcile),
        "all_stage5_reconciled": all(x["all_six_stage5_checks_pass"] for x in reconcile),
        "four_domains": DOMAINS,
        "unknown_indication_exact_terms": sorted(UNKNOWN_INDICATIONS),
        "exclusive_combinations_per_drug_database": 15,
        "nonexclusive_intersections_per_drug_database": 11,
        "all_exclusive_sums_equal_broad": all(
            r["broad_reports_n"] == r["exclusive_combinations_sum_n"] for r in broad_rows
        ),
        "all_intersections_reconstructed": True,
        "row_counts": {filename: len(rows) for filename, rows in outputs.items()},
        "canada_schema_url": SCHEMA_URL,
        "canada_schema_access_date": "2026-10-19",
        "dependency_script": str(
            Path(__file__).resolve().parent / "run_stage5_primary_analysis.py"
        ),
        "dependency_script_sha256": hashlib.sha256(
            (Path(__file__).resolve().parent / "run_stage5_primary_analysis.py").read_bytes()
        ).hexdigest(),
        "limitations": [
            "Recorded age/sex/indication does not validate patient-level confounders or treatment timing.",
            "Canada AGE_Y is numeric years; JADER age intervals and broad categories are retained without midpoint conversion.",
            "Missingness includes unknown separately; no imputation or report deletion.",
            "Nonexclusive domain contributions and intersections cannot be added; only the 15 exact patterns sum to the broad union.",
            "Report-level co-occurrence cannot establish clinical sequence, diagnosis confirmation, mechanism or causality.",
        ],
    }
    (args.out / "external_analysis_verification.json").write_text(
        json.dumps(final, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(final, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
