from __future__ import annotations
import csv, io, json, math, re, sqlite3, zipfile
from datetime import datetime
from pathlib import Path
import scipy
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
from statsmodels.stats.multitest import multipletests
import statsmodels

from package_paths import OUTPUT, FAERS_ZIPS, FAERS_DB, CANADA, JADER, STAGE3

TARGET = OUTPUT
OUT = TARGET / "03_确认后补充分析" / "阶段五_主要分析"
CAN = CANADA
STAGE3_FLOW = STAGE3 / "阶段三_报告级样本流程.csv"
DRUG_TERMS = {
    "alpelisib": ["ALPELISIB", "PIQRAY", "BYL719", "BYL-719", "アルペリシブ", "ピクレイ"],
    "inavolisib": ["INAVOLISIB", "ITOVEBI", "GDC0077", "GDC-0077", "イナボリシブ", "イトベビ"],
    "capivasertib": ["CAPIVASERTIB", "TRUQAP", "AZD5363", "AZD-5363", "カピバセルチブ", "トルカップ"],
    "everolimus": ["EVEROLIMUS", "AFINITOR", "RAD001", "RAD-001", "エベロリムス", "アフィニトール"],
    "temsirolimus": ["TEMSIROLIMUS", "TORISEL", "CCI779", "CCI-779", "テムシロリムス", "トーリセル"],
}
OUTCOMES = {
    "primary_union": ["Diabetic ketoacidosis", "Ketoacidosis", "Metabolic acidosis"],
    "PT_diabetic_ketoacidosis": ["Diabetic ketoacidosis"],
    "PT_ketoacidosis": ["Ketoacidosis"],
    "PT_metabolic_acidosis": ["Metabolic acidosis"],
    "hyperglycemia_domain": ["Hyperglycaemia", "Blood glucose increased", "Diabetes mellitus"],
    "dehydration_electrolyte_domain": [
        "Dehydration",
        "Hypovolaemia",
        "Hypokalaemia",
        "Hyponatraemia",
    ],
    "renal_injury_domain": ["Acute kidney injury", "Renal failure", "Renal impairment"],
    "severe_metabolic_renal_composite": [
        "Diabetic ketoacidosis",
        "Ketoacidosis",
        "Metabolic acidosis",
        "Dehydration",
        "Hypovolaemia",
        "Hypokalaemia",
        "Hyponatraemia",
        "Acute kidney injury",
        "Renal failure",
        "Renal impairment",
    ],
    "PT_acute_kidney_injury": ["Acute kidney injury"],
    "PT_renal_failure": ["Renal failure"],
    "PT_renal_impairment": ["Renal impairment"],
}
PRIMARY = "primary_union"
SECONDARY = [k for k in OUTCOMES if k != PRIMARY]
PT_EN = {
    "diabetic ketoacidosis",
    "ketoacidosis",
    "metabolic acidosis",
    "hyperglycaemia",
    "blood glucose increased",
    "diabetes mellitus",
    "dehydration",
    "hypovolaemia",
    "hypokalaemia",
    "hyponatraemia",
    "acute kidney injury",
    "renal failure",
    "renal impairment",
}
JP_TO_EN = {
    "糖尿病性ケトアシドーシス": "Diabetic ketoacidosis",
    "ケトアシドーシス": "Ketoacidosis",
    "代謝性アシドーシス": "Metabolic acidosis",
    "高血糖": "Hyperglycaemia",
    "血中ブドウ糖増加": "Blood glucose increased",
    "糖尿病": "Diabetes mellitus",
    "脱水": "Dehydration",
    "血液量減少症": "Hypovolaemia",
    "低カリウム血症": "Hypokalaemia",
    "低ナトリウム血症": "Hyponatraemia",
    "急性腎障害": "Acute kidney injury",
    "腎不全": "Renal failure",
    "腎機能障害": "Renal impairment",
}


def norm(s):
    return re.sub(r"[\W_]", "", (s or "").upper(), flags=re.UNICODE)


NORM_TERMS = {d: [norm(t) for t in ts] for d, ts in DRUG_TERMS.items()}


def match_drug(*fields):
    h = norm(" ".join(fields))
    for d, terms in NORM_TERMS.items():
        if any(t and t in h for t in terms):
            return d
    return None


def find_member(zf, prefix):
    for name in zf.namelist():
        base = Path(name).name.upper()
        if base.startswith(prefix.upper()) and base.endswith(".TXT"):
            return name
    return None


def iter_table(zf, member):
    with zf.open(member) as raw:
        txt = io.TextIOWrapper(raw, encoding="latin-1", errors="replace", newline="")
        for row in csv.DictReader(txt, delimiter="$"):
            yield {k.lower(): (v or "").strip() for k, v in row.items() if k is not None}


def write_csv(name, rows):
    if rows:
        with (OUT / name).open("w", encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for row in rows for k in row)))
            w.writeheader()
            w.writerows(rows)


def load_faers():
    zips = sorted(FAERS_ZIPS.glob("faers_ascii_*.zip"), key=lambda p: p.name.lower())
    zips = [p for p in zips if "2019q1" <= p.stem[-6:].lower() <= "2026q1"]
    if len(zips) != 29:
        raise RuntimeError(f"Expected 29 FAERS quarter zips; found {len(zips)}")
    conn = sqlite3.connect(FAERS_DB.resolve().as_uri() + "?mode=ro", uri=True)
    n, primaryids, caseids = conn.execute(
        "select count(*),count(distinct primaryid),count(distinct caseid) from demo_latest"
    ).fetchone()
    valid = {r[0] for r in conn.execute("select primaryid from demo_latest")}
    conn.close()
    if n != len(valid) or n != 10463182:
        raise RuntimeError(f"FAERS report count mismatch: {n}, loaded {len(valid)}")
    target = {d: set() for d in DRUG_TERMS}
    drug_rows = 0
    for zp in zips:
        with zipfile.ZipFile(zp) as zf:
            m = find_member(zf, "DRUG")
            if not m:
                continue
            for row in iter_table(zf, m):
                drug_rows += 1
                pid = row.get("primaryid", "")
                if pid in valid and row.get("role_cod", "").upper() == "PS":
                    d = match_drug(row.get("drugname", ""), row.get("prod_ai", ""))
                    if d:
                        target[d].add(pid)
    ptsets = {p: set() for p in PT_EN}
    reac_rows = 0
    for zp in zips:
        with zipfile.ZipFile(zp) as zf:
            m = find_member(zf, "REAC")
            if not m:
                continue
            for row in iter_table(zf, m):
                reac_rows += 1
                pid = row.get("primaryid", "")
                pt = row.get("pt", "").strip().lower()
                if pid in valid and pt in ptsets:
                    ptsets[pt].add(pid)
    events = {k: set().union(*(ptsets[p.lower()] for p in pts)) for k, pts in OUTCOMES.items()}
    meta = {
        "database": "FAERS",
        "background_reports_n": n,
        "primaryids_n": primaryids,
        "caseids_n": caseids,
        "drug_rows_scanned": drug_rows,
        "reaction_rows_scanned": reac_rows,
        "target_reports_by_drug": {k: len(v) for k, v in target.items()},
        "source_zips": [
            {
                "name": p.name,
                "bytes": p.stat().st_size,
                "modified": datetime.fromtimestamp(p.stat().st_mtime).isoformat(timespec="seconds"),
            }
            for p in zips
        ],
        "source_db": str(FAERS_DB),
    }
    return n, target, events, meta


def load_canada():
    reports = set()
    report_rows = 0
    with (CAN / "reports.txt").open(encoding="utf-8-sig", errors="replace", newline="") as f:
        for row in csv.reader(f, delimiter="$", quotechar='"'):
            if row:
                report_rows += 1
                reports.add(row[0].strip())
    target = {d: set() for d in DRUG_TERMS}
    drug_rows = 0
    with (CAN / "report_drug.txt").open(encoding="utf-8-sig", errors="replace", newline="") as f:
        for row in csv.reader(f, delimiter="$", quotechar='"'):
            drug_rows += 1
            if len(row) < 6 or row[4].strip().upper() != "SUSPECT":
                continue
            d = match_drug(row[3])
            if d and row[1].strip() in reports:
                target[d].add(row[1].strip())
    ptsets = {p: set() for p in PT_EN}
    reac_rows = 0
    with (CAN / "reactions.txt").open(encoding="utf-8-sig", errors="replace", newline="") as f:
        for row in csv.reader(f, delimiter="$", quotechar='"'):
            reac_rows += 1
            if len(row) < 6:
                continue
            rid = row[1].strip()
            pt = row[5].strip().lower()
            if rid in reports and pt in ptsets:
                ptsets[pt].add(rid)
    events = {k: set().union(*(ptsets[p.lower()] for p in pts)) for k, pts in OUTCOMES.items()}
    files = [CAN / "reports.txt", CAN / "report_drug.txt", CAN / "reactions.txt"]
    meta = {
        "database": "Canada Vigilance",
        "background_reports_n": len(reports),
        "report_rows_scanned": report_rows,
        "drug_rows_scanned": drug_rows,
        "reaction_rows_scanned": reac_rows,
        "target_reports_by_drug": {k: len(v) for k, v in target.items()},
        "source_files": [
            {
                "name": p.name,
                "bytes": p.stat().st_size,
                "modified": datetime.fromtimestamp(p.stat().st_mtime).isoformat(timespec="seconds"),
            }
            for p in files
        ],
    }
    return len(reports), target, events, meta


def jkey(rid, version):
    try:
        return f"{rid}|{int(version):02d}"
    except (TypeError, ValueError):
        return f"{rid}|{version}"


def load_jader():
    latest = {}
    demo_rows = 0
    with (JADER / "demo202605.csv").open(encoding="cp932", errors="replace", newline="") as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            if len(row) < 4:
                continue
            demo_rows += 1
            try:
                v = int(row[1])
            except ValueError:
                v = 0
            if row[0] not in latest or v > latest[row[0]]:
                latest[row[0]] = v
    reports = {jkey(rid, v) for rid, v in latest.items()}
    target = {d: set() for d in DRUG_TERMS}
    drug_rows = 0
    with (JADER / "drug202605.csv").open(encoding="cp932", errors="replace", newline="") as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            drug_rows += 1
            if len(row) < 6 or row[3].strip() != "被疑薬":
                continue
            key = jkey(row[0], row[1])
            if key not in reports:
                continue
            d = match_drug(row[4], row[5])
            if d:
                target[d].add(key)
    jp_to_pt = {jp: en.lower() for jp, en in JP_TO_EN.items()}
    ptsets = {p: set() for p in PT_EN}
    reac_rows = 0
    with (JADER / "reac202605.csv").open(encoding="cp932", errors="replace", newline="") as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            reac_rows += 1
            if len(row) < 4:
                continue
            key = jkey(row[0], row[1])
            pt = jp_to_pt.get(row[3].strip())
            if key in reports and pt in ptsets:
                ptsets[pt].add(key)
    events = {k: set().union(*(ptsets[p.lower()] for p in pts)) for k, pts in OUTCOMES.items()}
    files = [JADER / "demo202605.csv", JADER / "drug202605.csv", JADER / "reac202605.csv"]
    meta = {
        "database": "JADER",
        "background_reports_n": len(reports),
        "demo_rows_scanned": demo_rows,
        "superseded_report_versions_n": demo_rows - len(reports),
        "drug_rows_scanned": drug_rows,
        "reaction_rows_scanned": reac_rows,
        "target_reports_by_drug": {k: len(v) for k, v in target.items()},
        "source_files": [
            {
                "name": p.name,
                "bytes": p.stat().st_size,
                "modified": datetime.fromtimestamp(p.stat().st_mtime).isoformat(timespec="seconds"),
            }
            for p in files
        ],
    }
    return len(reports), target, events, meta


def fmt(x, digits=5):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return ""
    return f"{x:.{digits}g}"


def calc(db, drug, outcome, N, exposed, cases, category):
    a = len(exposed & cases)
    ab = len(exposed)
    ac = len(cases)
    b = ab - a
    c = ac - a
    d = N - a - b - c
    if min(a, b, c, d) < 0 or a + b + c + d != N:
        raise RuntimeError(f"invalid 2x2: {db}/{drug}/{outcome} {(a,b,c,d)} N={N}")
    estimable = ab > 0 and ac > 0 and N - ac > 0 and N - ab > 0
    if estimable:
        t = [[a, b], [c, d]]
        sample = (a * d) / (b * c) if b * c else (math.inf if a * d > 0 else math.nan)
        try:
            o = odds_ratio(t, kind="conditional")
            cond = float(o.statistic)
            ci = o.confidence_interval(confidence_level=0.95)
            lo = float(ci.low)
            hi = float(ci.high)
        except Exception:
            cond = lo = hi = math.nan
        p = float(fisher_exact(t, alternative="two-sided").pvalue)
    else:
        sample = cond = lo = hi = p = math.nan
    return {
        "database": db,
        "drug": drug,
        "outcome": outcome,
        "outcome_class": category,
        "a_target_case_reports": a,
        "b_target_noncase_reports": b,
        "c_background_case_reports": c,
        "d_background_noncase_reports": d,
        "target_reports_n": ab,
        "background_reports_n": N,
        "total_case_reports_n": ac,
        "target_case_report_share_pct": round(100 * a / ab, 4) if ab else "",
        "estimable": estimable,
        "ROR_cross_product": fmt(sample),
        "conditional_OR_MLE": fmt(cond),
        "exact_conditional_95CI_lower": fmt(lo),
        "exact_conditional_95CI_upper": fmt(hi),
        "fisher_two_sided_p": fmt(p),
        "_p": p,
    }


def adjust(rows, method, key):
    valid = [r for r in rows if r["estimable"]]
    if valid:
        vals = multipletests([r["_p"] for r in valid], method=method)[1]
        for r, v in zip(valid, vals):
            r[key] = fmt(float(v))


def clean(rows):
    return [{k: v for k, v in r.items() if not k.startswith("_")} for r in rows]


def read_locked():
    if not STAGE3_FLOW.is_file():
        return []
    with STAGE3_FLOW.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    loaded = []
    for loader in [load_faers, load_canada, load_jader]:
        print("Loading", loader.__name__, flush=True)
        item = loader()
        loaded.append(item)
        print(
            "Loaded",
            item[3]["database"],
            "background",
            item[0],
            "targets",
            item[3]["target_reports_by_drug"],
            flush=True,
        )
    data = {}
    qa = []
    meta_rows = []
    locked = read_locked()
    for N, target, events, meta in loaded:
        db = meta["database"]
        data[db] = (N, target, events)
        meta_rows.append(meta)
        qa.append(
            {
                "database": db,
                "background_reports_n": N,
                "all_database_primary_case_reports_n": len(events[PRIMARY]),
                "target_reports_by_drug_json": json.dumps(
                    meta["target_reports_by_drug"], ensure_ascii=False
                ),
                "event_reports_by_outcome_json": json.dumps(
                    {k: len(v) for k, v in events.items()}, ensure_ascii=False
                ),
            }
        )
        if not locked:
            if not all(len(target[d]) <= N and len(target[d] & events[PRIMARY]) <= len(target[d]) for d in DRUG_TERMS):
                raise ValueError(f"Invalid reconstructed report counts: {db}")
            continue
        for drug in DRUG_TERMS:
            e = next(r for r in locked if r["database"] == db and r["drug"] == drug)
            rn = len(target[drug])
            cn = len(target[drug] & events[PRIMARY])
            ok = (
                rn == int(e["target_suspect_reports_n"])
                and cn == int(e["primary_PT_union_reports_in_target_n"])
                and N == int(e["background_reports_n"])
            )
            qa.append(
                {
                    "database": db,
                    "drug": drug,
                    "rebuilt_target_reports_n": rn,
                    "stage3_target_reports_n": e["target_suspect_reports_n"],
                    "rebuilt_primary_case_reports_n": cn,
                    "stage3_primary_case_reports_n": e["primary_PT_union_reports_in_target_n"],
                    "background_reports_n": N,
                    "reconciles_stage3": ok,
                }
            )
            if not ok:
                raise RuntimeError(f"Stage 3 reconciliation failed: {db}/{drug}")
    main_rows = []
    sec_rows = []
    for db, (N, target, events) in data.items():
        for drug in DRUG_TERMS:
            main_rows.append(calc(db, drug, PRIMARY, N, target[drug], events[PRIMARY], "primary"))
            for out in SECONDARY:
                sec_rows.append(calc(db, drug, out, N, target[drug], events[out], "secondary"))
    adjust(main_rows, "holm", "holm_adjusted_p")
    adjust(sec_rows, "fdr_by", "BY_adjusted_q")
    for r in main_rows:
        lo = (
            float(r["exact_conditional_95CI_lower"])
            if r["exact_conditional_95CI_lower"]
            else math.nan
        )
        hp = float(r.get("holm_adjusted_p", "")) if r.get("holm_adjusted_p") else math.nan
        if not r["estimable"]:
            status = "not_estimable"
        elif r["a_target_case_reports"] >= 3 and lo > 1 and hp < 0.05:
            status = "passes_case_CI_Holm_rule"
        elif r["a_target_case_reports"] >= 3 and lo > 1:
            status = "nominal_screen_positive_not_Holm_significant"
        else:
            status = "does_not_pass_prespecified_rule"
        r["prespecified_primary_status"] = status
    write_csv("阶段五_队列核验.csv", qa)
    write_csv("阶段五_主要结局ROR_Fisher_Holm.csv", clean(main_rows))
    write_csv("阶段五_次要结局ROR_Fisher_BY.csv", clean(sec_rows))
    maprows = []
    for jp, en in JP_TO_EN.items():
        maprows.append(
            {
                "JADER_reaction_PT_exact": jp,
                "mapped_English_MedDRA_PT": en,
                "outcome_groups": ";".join(k for k, v in OUTCOMES.items() if en in v),
                "mapping_rule": "exact string match",
            }
        )
    write_csv("阶段五_JADER_PT映射.csv", maprows)
    write_csv(
        "阶段五_输入来源清单.csv",
        [
            {
                "source": "FAERS quarterly ASCII data",
                "path": str(FAERS_ZIPS),
                "files_n": 29,
                "window": "2019Q1-2026Q1",
                "database": "FAERS",
            },
            {
                "source": "FAERS latest-case DB (read-only)",
                "path": str(FAERS_DB),
                "database": "FAERS",
            },
            {
                "source": "Canada Vigilance",
                "path": str(CAN),
                "version": "cvponline_extract_20260228",
                "database": "Canada Vigilance",
            },
            {
                "source": "JADER",
                "path": str(JADER),
                "version": "demo202605/drug202605/reac202605",
                "database": "JADER",
            },
            {
                "source": "locked Stage 3 aggregate counts",
                "path": str(STAGE3_FLOW),
                "database": "all",
            },
        ],
    )
    nmain = sum(bool(r["estimable"]) for r in main_rows)
    nsec = sum(bool(r["estimable"]) for r in sec_rows)
    diagnostics = {
        "analysis_date": datetime.now().date().isoformat(),
        "python_packages": {"scipy": scipy.__version__, "statsmodels": statsmodels.__version__},
        "main_estimable_tests": nmain,
        "main_max_tests": len(main_rows),
        "secondary_BY_estimable_tests": nsec,
        "secondary_BY_max_tests": len(sec_rows),
        "interval": "SciPy conditional odds-ratio exact confidence interval; cross-product ROR separately reported",
        "p_value": "two-sided Fisher exact test",
        "multiplicity": "Holm for all estimable primary drug x database tests; BY across all estimable secondary tests",
        "BCPNN": "Not computed: no locally available independently verifiable implementation/reference test suite found; approximate IC/EBGM metrics not substituted.",
        "checks": {
            "prior_audit_reconciliation_performed": bool(locked),
            "stage3_reconciliation_all_passed": True if locked else None,
            "all_2x2_rows_sum_to_background": True,
        },
        "source_metadata": meta_rows,
    }
    (OUT / "阶段五_算法与核验记录.json").write_text(
        json.dumps(diagnostics, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    source_files = [
        *sorted(FAERS_ZIPS.glob("faers_ascii_*.zip")),
        FAERS_DB,
        CAN / "reports.txt",
        CAN / "report_drug.txt",
        CAN / "reactions.txt",
        JADER / "demo202605.csv",
        JADER / "drug202605.csv",
        JADER / "reac202605.csv",
        STAGE3_FLOW,
    ]
    write_csv(
        "阶段五_源文件指纹.csv",
        [
            {
                "path": str(p),
                "exists": p.exists(),
                "bytes": p.stat().st_size if p.exists() else "",
                "modified": datetime.fromtimestamp(p.stat().st_mtime).isoformat(timespec="seconds")
                if p.exists()
                else "",
            }
            for p in source_files
        ],
    )

    print(
        json.dumps(
            {
                "status": "completed",
                "out": str(OUT),
                "main_estimable": nmain,
                "secondary_estimable": nsec,
                "main": [
                    {
                        k: r.get(k, "")
                        for k in [
                            "database",
                            "drug",
                            "a_target_case_reports",
                            "ROR_cross_product",
                            "exact_conditional_95CI_lower",
                            "exact_conditional_95CI_upper",
                            "fisher_two_sided_p",
                            "holm_adjusted_p",
                            "prespecified_primary_status",
                        ]
                    }
                    for r in main_rows
                ],
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
