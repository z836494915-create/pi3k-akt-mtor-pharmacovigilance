from __future__ import annotations

import csv
import importlib.util
import io
import json
import re
import sqlite3
import zipfile
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import scipy
import statsmodels
from statsmodels.stats.multitest import multipletests


from package_paths import OUTPUT, FAERS_DB, FAERS_ZIPS, CANADA, STAGE3

ROOT = OUTPUT
STAGE5 = ROOT / "03_确认后补充分析" / "阶段五_主要分析"
OUT = ROOT / "03_确认后补充分析" / "阶段六_敏感性分析"
STAGE5_SCRIPT = Path(__file__).with_name("run_stage5_primary_analysis.py")
CAN = CANADA

ONCOLOGY_RE = re.compile(
    r"(cancer|carcinoma|neoplasm|tumou?r|sarcoma|melanoma|leukaemia|leukemia|lymphoma|myeloma|"
    r"blastoma|malignan|breast|prostate|lung|colorectal|renal cell|endometrial|ovarian|pancreatic|"
    r"glioblastoma|cholangiocarcinoma|hepatocellular)",
    re.I,
)


def load_stage5_module():
    spec = importlib.util.spec_from_file_location("stage5_primary", STAGE5_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import Stage 5 analysis definitions: {STAGE5_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_csv(path: Path, rows: list[dict]):
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def read_csv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def date_parse(raw: str):
    raw = (raw or "").strip()
    if not raw:
        return None
    for fmt in ("%Y%m%d", "%Y-%m-%d", "%d-%b-%Y", "%Y/%m/%d"):
        try:
            return datetime.strptime(raw[:10], fmt).date()
        except ValueError:
            pass
    return None


def read_faers_dates_and_ids():
    conn = sqlite3.connect(FAERS_DB.resolve().as_uri() + "?mode=ro", uri=True)
    cols = [r[1] for r in conn.execute("pragma table_info(demo_latest)")]
    by_upper = {c.upper(): c for c in cols}
    pid_col = by_upper.get("PRIMARYID")
    dt_col = by_upper.get("FDA_DT")
    if not pid_col or not dt_col:
        conn.close()
        raise RuntimeError(f"demo_latest does not have PRIMARYID/FDA_DT: {cols}")
    rows = conn.execute(f'SELECT "{pid_col}", "{dt_col}" FROM demo_latest')
    dates = {}
    for pid, raw in rows:
        if pid is not None:
            dates[str(pid)] = date_parse(str(raw or ""))
    conn.close()
    return dates


def find_member(zf, prefix):
    for name in zf.namelist():
        base = Path(name).name.upper()
        if base.startswith(prefix.upper()) and base.endswith(".TXT"):
            return name
    return None


def iter_table(zf, member):
    with zf.open(member) as raw:
        txt = io.TextIOWrapper(raw, encoding="latin-1", errors="replace", newline="")
        for row in csv.DictReader(txt, delimiter="$", quotechar='"'):
            yield {k.lower(): (v or "").strip() for k, v in row.items() if k is not None}


def load_faers_extended(s5):
    zips = sorted(FAERS_ZIPS.glob("faers_ascii_*.zip"), key=lambda p: p.name.lower())
    zips = [p for p in zips if "2019q1" <= p.stem[-6:].lower() <= "2026q1"]
    if len(zips) != 29:
        raise RuntimeError(f"Expected 29 FAERS quarterly files (2019Q1–2026Q1), found {len(zips)}")
    dates = read_faers_dates_and_ids()
    valid = set(dates)
    roles = {d: {"PS": set(), "SS": set()} for d in s5.DRUG_TERMS}
    ps_seq = defaultdict(set)
    drug_rows = 0
    for zp in zips:
        with zipfile.ZipFile(zp) as zf:
            member = find_member(zf, "DRUG")
            if not member:
                continue
            for row in iter_table(zf, member):
                drug_rows += 1
                pid = row.get("primaryid", "")
                if pid not in valid:
                    continue
                role = row.get("role_cod", "").upper()
                if role not in ("PS", "SS"):
                    continue
                drug = s5.match_drug(row.get("drugname", ""), row.get("prod_ai", ""))
                if not drug:
                    continue
                roles[drug][role].add(pid)
                if role == "PS":
                    seq = row.get("drug_seq", "").strip()
                    if seq:
                        ps_seq[(drug, pid)].add(seq)

    targetseqs_by_pid = defaultdict(dict)
    for (drug, pid), seqs in ps_seq.items():
        targetseqs_by_pid[pid][drug] = seqs
    linked_any = {d: set() for d in s5.DRUG_TERMS}
    linked_oncology = {d: set() for d in s5.DRUG_TERMS}
    oncology = set()
    indi_rows = 0
    for zp in zips:
        with zipfile.ZipFile(zp) as zf:
            member = find_member(zf, "INDI")
            if not member:
                continue
            for row in iter_table(zf, member):
                indi_rows += 1
                pid = row.get("primaryid", "")
                if pid not in valid:
                    continue
                indi = row.get("indi_pt", "")
                is_oncology = bool(ONCOLOGY_RE.search(indi))
                if is_oncology:
                    oncology.add(pid)
                seq = row.get("indi_drug_seq", "").strip()
                if not seq:
                    continue
                for drug, seqs in targetseqs_by_pid.get(pid, {}).items():
                    if seq in seqs:
                        linked_any[drug].add(pid)
                        if is_oncology:
                            linked_oncology[drug].add(pid)

    ps = {d: roles[d]["PS"] for d in s5.DRUG_TERMS}
    ps_union = set().union(*ps.values()) if ps else set()
    members = defaultdict(list)
    for drug, ids in ps.items():
        for pid in ids:
            members[pid].append(drug)
    ptsets = {p.lower(): set() for p in s5.PT_EN}
    target_ptsets = {d: defaultdict(set) for d in s5.DRUG_TERMS}
    reac_rows = 0
    for zp in zips:
        with zipfile.ZipFile(zp) as zf:
            member = find_member(zf, "REAC")
            if not member:
                continue
            for row in iter_table(zf, member):
                reac_rows += 1
                pid = row.get("primaryid", "")
                if pid not in valid:
                    continue
                pt = row.get("pt", "").strip().lower()
                if pt in ptsets:
                    ptsets[pt].add(pid)
                if pid in ps_union and pt:
                    for drug in members[pid]:
                        target_ptsets[drug][pt].add(pid)
    events = {k: set().union(*(ptsets[p.lower()] for p in pts)) for k, pts in s5.OUTCOMES.items()}
    meta = {
        "background_reports_n": len(valid),
        "valid_FDA_DT_reports_n": sum(x is not None for x in dates.values()),
        "dates": dates,
        "roles": roles,
        "PS": ps,
        "PS_SS": {d: roles[d]["PS"] | roles[d]["SS"] for d in s5.DRUG_TERMS},
        "SS_only": {d: roles[d]["SS"] - roles[d]["PS"] for d in s5.DRUG_TERMS},
        "role_overlap": {d: roles[d]["SS"] & roles[d]["PS"] for d in s5.DRUG_TERMS},
        "oncology": oncology,
        "linked_any": linked_any,
        "linked_oncology": linked_oncology,
        "target_ptsets": target_ptsets,
        "source_zips": zips,
        "counts": {"drug_rows": drug_rows, "indi_rows": indi_rows, "reaction_rows": reac_rows},
    }
    return valid, events, meta


def run_calc(s5, sensitivity, db, drug, universe, exposed, cases):
    exposed = exposed & universe
    cases = cases & universe
    if not exposed:
        return {
            "sensitivity": sensitivity,
            "sensitivity_family": sensitivity,
            "database": db,
            "drug": drug,
            "outcome": "primary_union",
            "estimable": False,
            "reason_not_estimable": "no exposed reports in restricted cohort",
            "cohort_reports_n": len(universe),
            "target_reports_n": 0,
        }
    r = s5.calc(db, drug, "primary_union", len(universe), exposed, cases, "sensitivity")
    r["sensitivity"] = sensitivity
    r["sensitivity_family"] = sensitivity
    r["cohort_reports_n"] = len(universe)
    r["reason_not_estimable"] = "" if r["estimable"] else "zero target/event or comparator cell"
    return r


def adjust_sensitivities(rows):
    groups = defaultdict(list)
    for r in rows:
        if r.get("estimable") and r.get("_p") is not None:
            groups[r["sensitivity_family"]].append(r)
    for family, group in groups.items():
        q = multipletests([r["_p"] for r in group], method="fdr_by")[1]
        for row, value in zip(group, q):
            row["BY_adjusted_q_within_family"] = f"{float(value):.7g}"
    for row in rows:
        row.setdefault("BY_adjusted_q_within_family", "")


def clean(rows):
    return [{k: v for k, v in r.items() if not k.startswith("_")} for r in rows]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    s5 = load_stage5_module()
    print(
        "Scanning FAERS roles, indications, reactions and FDA report dates (read-only)...",
        flush=True,
    )
    faers_universe, faers_events, fm = load_faers_extended(s5)
    print(
        f"FAERS scan complete: {len(faers_universe):,} reports; targets { {d: len(fm['PS'][d]) for d in s5.DRUG_TERMS} }",
        flush=True,
    )
    print("Loading Stage 5 Canada Vigilance and JADER report-level sets...", flush=True)
    can_n, can_target, can_events, can_meta = s5.load_canada()
    jad_n, jad_target, jad_events, jad_meta = s5.load_jader()
    data = {
        "FAERS": (faers_universe, {d: fm["PS"][d] for d in s5.DRUG_TERMS}, faers_events),
        "Canada Vigilance": (set(), can_target, can_events),
        "JADER": (set(), jad_target, jad_events),
    }

    # Canada report-type field is the English report type at zero-based column 7.
    canada_spontaneous = set()
    with (CAN / "reports.txt").open(encoding="utf-8-sig", errors="replace", newline="") as f:
        for row in csv.reader(f, delimiter="$", quotechar='"'):
            if len(row) > 7 and row[7].strip().casefold() == "spontaneous":
                canada_spontaneous.add(row[0].strip())

    sens_rows = []
    # FAERS role definition: broaden PS to PS+SS.
    for drug in s5.DRUG_TERMS:
        sens_rows.append(
            run_calc(
                s5,
                "FAERS_PS+SS",
                "FAERS",
                drug,
                faers_universe,
                fm["PS_SS"][drug],
                faers_events["primary_union"],
            )
        )
        # Keep PS as exposed, but remove SS-only reports of the same target drug from comparators.
        u = faers_universe - fm["SS_only"][drug]
        sens_rows.append(
            run_calc(
                s5,
                "FAERS_PS_vs_no_target_SS_only",
                "FAERS",
                drug,
                u,
                fm["PS"][drug],
                faers_events["primary_union"],
            )
        )

    # Restrict the report universe to Canada spontaneous reports.
    for drug in s5.DRUG_TERMS:
        sens_rows.append(
            run_calc(
                s5,
                "Canada_spontaneous_only",
                "Canada Vigilance",
                drug,
                canada_spontaneous,
                can_target[drug],
                can_events["primary_union"],
            )
        )

    # Cancer-report restriction uses any oncology-like INDI PT for the background;
    # target exposure additionally requires exact DRUG.drug_seq = INDI.indi_drug_seq linkage.
    for drug in s5.DRUG_TERMS:
        exp = fm["linked_oncology"][drug] & fm["PS"][drug]
        u = fm["oncology"] - (fm["PS_SS"][drug] - exp)
        sens_rows.append(
            run_calc(
                s5,
                "FAERS_oncology_INDI_exact_target_sequence",
                "FAERS",
                drug,
                u,
                exp,
                faers_events["primary_union"],
            )
        )

    # FDA_DT is report-received date; report quarters are counts only, annual strata are sensitivity tables.
    years = defaultdict(set)
    quarters = defaultdict(set)
    for pid, day in fm["dates"].items():
        if day:
            years[day.year].add(pid)
            quarters[f"{day.year}Q{(day.month - 1) // 3 + 1}"].add(pid)
    annual_rows = []
    quarterly_counts = []
    for quarter, universe in sorted(quarters.items()):
        cases = faers_events["primary_union"] & universe
        for drug in s5.DRUG_TERMS:
            exposed = fm["PS"][drug] & universe
            quarterly_counts.append(
                {
                    "database": "FAERS",
                    "report_quarter_by_FDA_DT": quarter,
                    "drug": drug,
                    "quarter_reports_n": len(universe),
                    "PS_target_reports_n": len(exposed),
                    "primary_union_reports_in_quarter_n": len(cases),
                    "primary_union_in_PS_target_reports_n": len(exposed & cases),
                }
            )
    for year, universe in sorted(years.items()):
        for drug in s5.DRUG_TERMS:
            row = run_calc(
                s5,
                "FAERS_annual_by_FDA_DT",
                "FAERS",
                drug,
                universe,
                fm["PS"][drug],
                faers_events["primary_union"],
            )
            row["report_year_by_FDA_DT"] = year
            annual_rows.append(row)
    sens_rows.extend(annual_rows)

    adjust_sensitivities(sens_rows)

    # Report-level decomposition and co-report counts in all three databases.
    decom_rows = []
    coreport_rows = []
    datasets = {
        "FAERS": (faers_universe, fm["PS"], faers_events),
        "Canada Vigilance": (set(range(can_n)), can_target, can_events),
        "JADER": (set(range(jad_n)), jad_target, jad_events),
    }
    for db, (universe_marker, targets, events) in datasets.items():
        for drug in s5.DRUG_TERMS:
            exposed = targets[drug]
            primary = events["primary_union"]
            dehy = events["dehydration_electrolyte_domain"]
            renal = events["renal_injury_domain"]
            broad = events["severe_metabolic_renal_composite"]
            hyper = events["hyperglycemia_domain"]
            pids = exposed & primary
            decom_rows.append(
                {
                    "database": db,
                    "drug": drug,
                    "PS_or_suspect_target_reports_n": len(exposed),
                    "PT_DKA_reports_n": len(exposed & events["PT_diabetic_ketoacidosis"]),
                    "PT_Ketoacidosis_reports_n": len(exposed & events["PT_ketoacidosis"]),
                    "PT_Metabolic_acidosis_reports_n": len(
                        exposed & events["PT_metabolic_acidosis"]
                    ),
                    "primary_union_reports_n": len(pids),
                    "dehydration_electrolyte_domain_reports_n": len(exposed & dehy),
                    "renal_injury_domain_reports_n": len(exposed & renal),
                    "primary_and_dehydration_electrolyte_overlap_n": len(exposed & primary & dehy),
                    "primary_and_renal_overlap_n": len(exposed & primary & renal),
                    "broad_composite_reports_n": len(exposed & broad),
                    "broad_composite_reports_outside_primary_n": len(exposed & (broad - primary)),
                }
            )
            coreport_rows.append(
                {
                    "database": db,
                    "drug": drug,
                    "PS_or_suspect_target_reports_n": len(exposed),
                    "primary_union_n": len(pids),
                    "hyperglycemia_domain_n": len(exposed & hyper),
                    "primary_hyperglycemia_overlap_n": len(exposed & primary & hyper),
                    "primary_cases_with_hyperglycemia_pct": round(
                        100 * len(exposed & primary & hyper) / len(pids), 2
                    )
                    if pids
                    else "",
                    "dehydration_electrolyte_domain_n": len(exposed & dehy),
                    "primary_dehydration_electrolyte_overlap_n": len(exposed & primary & dehy),
                    "renal_injury_domain_n": len(exposed & renal),
                    "primary_renal_overlap_n": len(exposed & primary & renal),
                }
            )

    # FAERS top five non-primary preferred terms by unique PS report count, then overlap with primary outcome.
    toppt_rows = []
    primary_pts = {x.lower() for x in s5.OUTCOMES["primary_union"]}
    for drug in s5.DRUG_TERMS:
        ranked = [
            (pt, ids) for pt, ids in fm["target_ptsets"][drug].items() if pt not in primary_pts
        ]
        ranked.sort(key=lambda x: (-len(x[1]), x[0]))
        for rank, (pt, ids) in enumerate(ranked[:5], start=1):
            toppt_rows.append(
                {
                    "database": "FAERS",
                    "drug": drug,
                    "rank_by_unique_target_report_count": rank,
                    "reaction_PT": pt,
                    "PS_target_reports_with_PT_n": len(ids),
                    "reports_with_PT_and_primary_union_n": len(ids & faers_events["primary_union"]),
                    "PT_reports_outside_primary_union_n": len(ids - faers_events["primary_union"]),
                    "interpretation": "co-report diagnostic only; no event removal or formal masking correction",
                }
            )

    # Canada spontaneous field/data linkage audit.
    can_type_counts = Counter()
    with (CAN / "reports.txt").open(encoding="utf-8-sig", errors="replace", newline="") as f:
        for row in csv.reader(f, delimiter="$", quotechar='"'):
            if row:
                can_type_counts[
                    (row[7].strip() if len(row) > 7 and row[7].strip() else "[blank]")
                ] += 1
    canada_rows = [
        {"report_type": k, "reports_n": v, "sensitivity_included": k.casefold() == "spontaneous"}
        for k, v in sorted(can_type_counts.items(), key=lambda x: (-x[1], x[0]))
    ]

    # Sequence-link audit: missing/ambiguous links are not classified as negative indication.
    indication_rows = []
    for drug in s5.DRUG_TERMS:
        target = fm["PS"][drug]
        linked = fm["linked_any"][drug]
        oncolink = fm["linked_oncology"][drug]
        onco_reports = fm["oncology"]
        indication_rows.append(
            {
                "database": "FAERS",
                "drug": drug,
                "PS_target_reports_n": len(target),
                "PS_target_reports_with_exact_drug_seq_INDI_link_n": len(target & linked),
                "PS_target_reports_with_exact_linked_oncology_INDI_n": len(target & oncolink),
                "all_oncology_INDI_background_reports_n": len(onco_reports),
                "PS_target_reports_without_exact_link_n": len(target - linked),
                "interpretation": "no exact linked oncology indication is not evidence of no cancer indication; missing/uncertain links remain unclassified",
            }
        )

    # Persist Stage 3 TTO/date and duplicate limitations already audited; do not infer/merge records.
    avail_path = STAGE3 / "阶段三_FAERS人口学与结局资料可用性.csv"
    availability = read_csv(avail_path) if avail_path.exists() else []
    tto_rows = [
        r
        for r in availability
        if any(k in " ".join(r.keys()).lower() for k in ("tto", "therapy", "event date", "caseid"))
    ]

    qa = []
    stage5_main = read_csv(STAGE5 / "阶段五_主要结局ROR_Fisher_Holm.csv")
    for drug in s5.DRUG_TERMS:
        st5 = next(r for r in stage5_main if r["database"] == "FAERS" and r["drug"] == drug)
        n_ps = len(fm["PS"][drug])
        c_ps = len(fm["PS"][drug] & faers_events["primary_union"])
        qa.append(
            {
                "check": "FAERS Stage 5 reconciliation",
                "drug": drug,
                "Stage5_target_reports_n": st5["target_reports_n"],
                "Stage6_PS_target_reports_n": n_ps,
                "Stage5_primary_cases_n": st5["a_target_case_reports"],
                "Stage6_primary_cases_n": c_ps,
                "passed": int(st5["target_reports_n"]) == n_ps
                and int(st5["a_target_case_reports"]) == c_ps,
            }
        )
    if not all(r["passed"] for r in qa):
        raise RuntimeError("FAERS target or primary-case counts do not reconcile with Stage 5")
    qa.append(
        {
            "check": "Canada spontaneous report-type count",
            "drug": "all",
            "Stage5_target_reports_n": "",
            "Stage6_PS_target_reports_n": len(canada_spontaneous),
            "Stage5_primary_cases_n": "",
            "Stage6_primary_cases_n": len(canada_spontaneous),
            "passed": len(canada_spontaneous) == 681308,
        }
    )

    write_csv(OUT / "阶段六_敏感性ROR结果.csv", clean(sens_rows))
    write_csv(OUT / "阶段六_FAERS季度报告计数.csv", quarterly_counts)
    write_csv(OUT / "阶段六_复合结局分解.csv", decom_rows)
    write_csv(OUT / "阶段六_共报矩阵.csv", coreport_rows)
    write_csv(OUT / "阶段六_FAERS高频非主要PT共报诊断.csv", toppt_rows)
    write_csv(OUT / "阶段六_FAERS适应证序号链接核验.csv", indication_rows)
    write_csv(OUT / "阶段六_Canada报告类型核验.csv", canada_rows)
    write_csv(OUT / "阶段六_TTO与重复报告可行性.csv", tto_rows)
    write_csv(OUT / "阶段六_核验.csv", qa)

    source_paths = [
        FAERS_DB,
        *fm["source_zips"],
        CAN / "reports.txt",
        CAN / "report_drug.txt",
        CAN / "reactions.txt",
        STAGE5_SCRIPT,
        STAGE5 / "阶段五_主要结局ROR_Fisher_Holm.csv",
        avail_path,
    ]
    source_rows = []
    for p in source_paths:
        source_rows.append(
            {
                "path": str(p),
                "exists": p.exists(),
                "bytes": p.stat().st_size if p.exists() else "",
                "modified": datetime.fromtimestamp(p.stat().st_mtime).isoformat(timespec="seconds")
                if p.exists()
                else "",
            }
        )
    write_csv(OUT / "阶段六_输入来源指纹.csv", source_rows)

    # Remove bulky in-memory date dictionary before serializing diagnostics.
    sens_by_family = Counter(r["sensitivity_family"] for r in sens_rows)
    diagnostics = {
        "analysis_date": datetime.now().date().isoformat(),
        "python_packages": {"scipy": scipy.__version__, "statsmodels": statsmodels.__version__},
        "read_only_sources": True,
        "FAERS": {
            "report_n": len(faers_universe),
            "valid_FDA_DT_n": fm["counts"].get(
                "valid_FDA_DT_n", sum(x is not None for x in fm["dates"].values())
            ),
            "row_scan_counts": fm["counts"],
            "target_role_counts": {
                d: {
                    "PS": len(fm["roles"][d]["PS"]),
                    "SS": len(fm["roles"][d]["SS"]),
                    "SS_only": len(fm["SS_only"][d]),
                    "PS_and_SS": len(fm["role_overlap"][d]),
                }
                for d in s5.DRUG_TERMS
            },
            "oncology_INDI_background_n": len(fm["oncology"]),
            "annual_periods": sorted(years),
            "quarter_periods": sorted(quarters),
        },
        "Canada": {
            "all_reports_n": can_n,
            "spontaneous_reports_n": len(canada_spontaneous),
            "target_reports_by_drug": {d: len(can_target[d]) for d in s5.DRUG_TERMS},
        },
        "JADER": {
            "all_reports_n": jad_n,
            "target_reports_by_drug": {d: len(jad_target[d]) for d in s5.DRUG_TERMS},
        },
        "sensitivity_tests_by_family": dict(sens_by_family),
        "multiplicity": "Benjamini-Yekutieli q within each prespecified sensitivity family; descriptive co-report tables have no p-values.",
        "masking": "Report-level co-report diagnostics only. No cases/events were removed and no formal masking correction was estimated.",
        "indication_restriction": "FAERS oncology-like INDI background, with exact target DRUG.drug_seq to INDI.indi_drug_seq for exposed target reports. Regex classification is exploratory; absent/mismatched links are unclassified.",
        "time": "FAERS FDA_DT (report receipt date), not event date; quarterly counts and annual stratum estimates. 2026 is partial.",
        "duplicate_reports": "No automated cross-CASEID merging: no stable patient identifier or validated duplicate rule available.",
        "tto": "Retained Stage 3 date/THER completeness audit only; TTO is not used as an outcome or adjusted endpoint.",
        "checks": {
            "faers_stage5_target_and_case_counts_all_passed": True,
            "canada_spontaneous_count_matches_stage3": len(canada_spontaneous) == 681308,
            "all_sensitivity_tests_recomputed_from_report_sets": True,
        },
        "source_metadata": [can_meta, jad_meta],
    }
    (OUT / "阶段六_算法与核验记录.json").write_text(
        json.dumps(diagnostics, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print(f"Stage 6 outputs written to: {OUT}", flush=True)


if __name__ == "__main__":
    main()
