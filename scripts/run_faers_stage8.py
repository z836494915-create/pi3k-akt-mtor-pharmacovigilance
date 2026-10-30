from __future__ import annotations
import csv, hashlib, importlib.util, itertools, json, math, pickle, re, sqlite3, sys, time, zipfile
from collections import Counter, defaultdict
from datetime import datetime
from functools import lru_cache
from pathlib import Path
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
from statsmodels.stats.multitest import multipletests

from package_paths import OUTPUT

ROOT = OUTPUT
OUT = OUTPUT / "03_确认后补充分析" / "阶段八_核验后修复" / "FAERS"
S5 = Path(__file__).resolve().parent / "run_stage5_primary_analysis.py"
spec = importlib.util.spec_from_file_location("s5", S5)
s5 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s5)
DRUGS = list(s5.DRUG_TERMS)
PTS = [
    "diabetic ketoacidosis",
    "ketoacidosis",
    "metabolic acidosis",
    "dehydration",
    "hypovolaemia",
    "hypokalaemia",
    "hyponatraemia",
    "acute kidney injury",
    "renal failure",
    "renal impairment",
    "hyperglycaemia",
    "blood glucose increased",
    "diabetes mellitus",
]
PTBITS = {x: 1 << i for i, x in enumerate(PTS)}
ONCO = re.compile(
    r"cancer|carcinoma|neoplasm|tumou?r|sarcoma|melanoma|leukaemia|leukemia|lymphoma|myeloma|blastoma|malignan|breast|prostate|lung|colorectal|renal cell|endometrial|ovarian|pancreatic|glioblastoma|cholangiocarcinoma|hepatocellular",
    re.I,
)
BREAST = re.compile(
    r"breast.*(?:cancer|carcinoma|neoplasm|malignan|tumou?r)|(?:cancer|carcinoma|neoplasm|malignan|tumou?r).*breast",
    re.I,
)
RENAL = re.compile(
    r"renal cell carcinoma|renal carcinoma|kidney carcinoma|kidney cancer|renal cancer", re.I
)
AGE_FACTORS = {"YR": 1.0, "MON": 1 / 12, "WK": 1 / 52.1775, "DY": 1 / 365.25, "DEC": 10.0}
CATS = {
    "glucocorticoid_ingredient": "dexamethasone prednisone prednisolone hydrocortisone methylprednisolone budesonide betamethasone triamcinolone cortisone deflazacort",
    "endocrine_anticancer": "fulvestrant letrozole anastrozole exemestane tamoxifen toremifene leuprolide leuprorelin goserelin triptorelin degarelix relugolix bicalutamide enzalutamide apalutamide darolutamide abiraterone megestrol",
    "glucose_lowering": "metformin glipizide glimepiride glyburide glibenclamide gliclazide repaglinide nateglinide pioglitazone rosiglitazone acarbose miglitol sitagliptin saxagliptin linagliptin alogliptin vildagliptin exenatide liraglutide dulaglutide semaglutide lixisenatide tirzepatide dapagliflozin empagliflozin canagliflozin ertugliflozin sotagliflozin",
    "selected_other_anticancer": "palbociclib ribociclib abemaciclib paclitaxel docetaxel doxorubicin epirubicin cyclophosphamide carboplatin cisplatin fluorouracil capecitabine gemcitabine vinorelbine eribulin trastuzumab pertuzumab bevacizumab pembrolizumab nivolumab ipilimumab atezolizumab avelumab durvalumab axitinib sunitinib pazopanib cabozantinib lenvatinib sorafenib regorafenib tivozanib olaparib talazoparib rucaparib niraparib lapatinib neratinib tucatinib imatinib dasatinib nilotinib bosutinib ponatinib gefitinib erlotinib afatinib osimertinib crizotinib alectinib brigatinib lorlatinib rituximab bendamustine bortezomib carfilzomib ixazomib lenalidomide pomalidomide thalidomide azacitidine decitabine cytarabine etoposide bleomycin methotrexate vincristine vinblastine oxaliplatin irinotecan temozolomide",
}
CAT_TERMS = {k: set(v.upper().split()) for k, v in CATS.items()}
INDICATION_UNKNOWN_TERMS = {
    "product used for unknown indication",
    "drug use for unknown indication",
    "unknown indication",
    "indication unknown",
    "unknown",
    "unspecified indication",
    "no indication",
    "not reported",
    "not provided",
    "not stated",
    "not known",
    "n/a",
}
PLAN = {
    "created_at": datetime.now().isoformat(timespec="seconds"),
    "analysis_type": "Post-review exploratory report-level case-noncase analyses; original primary analysis unchanged",
    "unit": "Latest eligible FAERS CASEID version as defined by existing audited demo_latest; not unique patients",
    "sources": "Same 29 quarterly releases 2019Q1–2026Q1 and 10,463,182 latest reports",
    "primary_PT_union": PTS[:3],
    "sole_suspect_rule": "Exactly one distinct nonmissing normalized PROD_AI component among ALL PS/SS records in BOTH target and comparator; FDA backslash and semicolon component separators split. Ambiguous /+| or AND formulations, blank fields and explicit UNKNOWN/UNSPECIFIED INGREDIENT placeholders excluded. Target aliases collapsed to target generic; other salt/hydrate/spelling variants conservatively remain distinct, so some same-parent products can be excluded. This is a restriction to unambiguous single encoded suspect ingredient, not evidence of monotherapy.",
    "indication_rule": "Target exposure requires exact target PS DRUG_SEQ–INDI_DRUG_SEQ linkage to category. Background requires linkage of any PS/SS drug sequence to same category. Target PS reports lacking target-specific category are removed, not reassigned to comparator. Breast text / renal cancer text dictionaries are descriptive and not cancer stage validation.",
    "indication_regexes": {"breast_cancer": BREAST.pattern, "renal_cancer": RENAL.pattern},
    "demographics": "Among target PS primary-outcome reports: age 0<=years<=120 and supported units; groups <18,18-44,45-64,65-74,>=75,unusable. Sex F/M/unknown retained explicitly. Linked indications reported non-exclusively; missing denominator is distinct reports, not indication rows.",
    "temporal": "FDA_DT receipt dates, all observed complete dates: yearly and quarterly four-cells. Show ALL estimable strata, not only significant strata. Fixed comparable full calendar window 2025-01-01 through 2025-12-31; window 2024-01-01 through 2025-12-31 secondary. Both are receipt-date windows, not equal time since market entry. Partial2026Q1 separately labeled.",
    "duplicate_rule": "Stringent cross-CASEID candidate groups require identical valid full EVENT_DT and INIT_FDA_DT, exact normalized reported age+unit, F/M sex, nonblank occurrence country (reporter country fallback), available recorded active-ingredient+role sets and PT sets; missing/ambiguous active ingredients disqualify a candidate. Recorded sets do not prove clinical or coding completeness; blank PTs are not counted and role completeness has not been validated. Keep smallest CASEID only for sensitivity. These are candidates, not proven duplicate patients. Cohort restriction applied symmetrically.",
    "masking_rule": "Before new results, statistical review replaced outcome-dependent comparator deletion with symmetric report restrictions: remove ALL hyperglycemia-domain reports; alternatively remove ALL diarrhea-PT reports; for everolimus/temsirolimus also remove ALL reports with PS/SS alpelisib, inavolisib or capivasertib, including co-reports with the target mTOR drug. Same restriction to target and comparator. These are co-report/background sensitivity analyses, not validated masking correction.",
    "broad_decomposition": {
        "primary": PTS[:3],
        "dehydration": PTS[3:5],
        "electrolyte": PTS[5:7],
        "renal": PTS[7:10],
    },
    "multiplicity": "One joint BY family for new non-temporal sensitivity tests; yearly and quarterly tests form two separately labeled exploratory BY families. No significance-driven selection. Descriptive intersections have no hypothesis tests.",
    "drug_classification": "Curated active-ingredient dictionary, non-exhaustive. Any-role co-record and PS/SS co-record are separately counted. Insulin* ingredients additionally classified glucose-lowering. Routes/timing unavailable; absence of dictionary match is not absence of use.",
    "audit_gate": "Stop rather than change definitions if audit/baseline reconciliation fails.",
}


def csvout(name, rows):
    if not rows:
        return
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        w.writeheader()
        w.writerows(rows)


def jwrite(name, data):
    (OUT / name).write_text(
        json.dumps(data, ensure_ascii=False, indent=2, default=str), encoding="utf-8"
    )


@lru_cache(maxsize=100000)
def validdate(x):
    x = str(x or "").strip()
    if not re.fullmatch(r"\d{8}", x):
        return 0
    try:
        datetime.strptime(x, "%Y%m%d")
        return int(x)
    except ValueError:
        return 0


def ageyears(a, u):
    try:
        x = float(a) * AGE_FACTORS.get(str(u).upper(), math.nan)
        return x if math.isfinite(x) and 0 <= x <= 120 else None
    except (TypeError, ValueError):
        return None


def norm(x):
    return re.sub(r"[^A-Z0-9]", "", str(x).upper())


def hash64(x):
    return int.from_bytes(hashlib.blake2b(x.encode(), digest_size=8).digest(), "little") or 1


def ingredients(ai):
    ai = ai.upper().strip()
    if not ai or re.search(r"[/+|]|\bAND\b", ai):
        return [], False
    comps = []
    for part in re.split(r"[;\\]", ai):
        part = part.strip()
        if not part:
            return [], False
        if not norm(part):
            return [], False
        if norm(part) in {
            "UNKNOWN",
            "UNKNOWNDRUG",
            "UNKNOWNINGREDIENT",
            "UNSPECIFIED",
            "UNSPECIFIEDDRUG",
            "UNSPECIFIEDINGREDIENT",
            "NA",
            "NONE",
            "NULL",
            "UNK",
            "NOTREPORTED",
        }:
            return [], False
        d = s5.match_drug(part)
        comps.append(d.upper() if d else norm(part))
    return sorted(set(comps)), bool(comps)


def classtags(ai):
    # Token boundary matching prevents e.g. a drug substring inside a different ingredient.
    tokens = set(re.findall(r"[A-Z]+", ai.upper()))
    return sum(
        (1 << i)
        for i, (k, vs) in enumerate(CAT_TERMS.items())
        if tokens & vs
        or (
            k == "glucose_lowering"
            and any(t.startswith("INSULIN") for t in tokens)
            and not re.search(r"INSULIN[- ]LIKE", ai.upper())
        )
    )


def known_indication(text):
    value = re.sub(r"\s+", " ", text.strip().lower())
    return (
        bool(value)
        and value not in INDICATION_UNKNOWN_TERMS
        and not re.fullmatch(
            r"(?:drug|product|medication)(?: use| used)? for (?:an? )?(?:unknown|unspecified) indication",
            value,
        )
    )


def rows(zf, prefix, counts):
    m = s5.find_member(zf, prefix)
    if not m:
        raise ValueError(f"Missing {prefix} in {zf.filename}")
    with zf.open(m) as f:
        header = f.readline().decode("latin-1").rstrip("\r\n").split("$")
        h = {x.lower(): i for i, x in enumerate(header)}
        for ln in f:
            counts[prefix + "_raw_rows"] += 1
            cells = ln.rstrip(b"\r\n").split(b"$")
            if len(cells) < len(header):
                counts[prefix + "_short_rows"] += 1
                continue
            try:
                pid = int(cells[h["primaryid"]].strip(b'" '))
            except ValueError:
                counts[prefix + "_invalid_pid"] += 1
                continue
            yield pid, cells, h


def val(c, h, k):
    return c[h[k]].strip(b'" ').decode("latin-1") if k in h else ""


def extract():
    OUT.mkdir(parents=True, exist_ok=True)
    if (OUT / "analysis_plan_locked.json").exists():
        previous = json.loads((OUT / "analysis_plan_locked.json").read_text(encoding="utf-8"))
        changed = [k for k in PLAN if k != "created_at" and previous.get(k) != PLAN[k]]
        if changed:
            stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            jwrite("analysis_plan_previous_before_results_" + stamp + ".json", previous)
            jwrite(
                "analysis_plan_change_before_results_" + stamp + ".json",
                {
                    "time": datetime.now().isoformat(),
                    "reason": "Before extraction complete and before any new effect estimate: independent statistical review requires symmetric co-report/background restrictions; raw PROD_AI audit identified backslash ingredient separators and same-parent salt/hydrate aliases, requiring conservative component parsing. No new effect results were viewed.",
                    "changed_fields": changed,
                },
            )
        else:
            PLAN["created_at"] = previous.get("created_at", PLAN["created_at"])
    jwrite("analysis_plan_locked.json", PLAN)
    jwrite(
        "drug_class_dictionary.json",
        {
            "classes": CATS,
            "additional_rule": "INSULIN-prefixed word => glucose_lowering",
            "scope": "non-exhaustive active-ingredient co-record classification; routes and timing unknown",
        },
    )
    cn = sqlite3.connect(s5.FAERS_DB.resolve().as_uri() + "?mode=ro", uri=True)
    n, npid, ncase = cn.execute(
        "select count(*),count(distinct primaryid),count(distinct caseid) from demo_latest"
    ).fetchone()
    if (n, npid, ncase) != (10463182, 10463182, 10463182):
        raise RuntimeError("DEMO baseline uniqueness mismatch")
    print("Loading compact report index", flush=True)
    pids = np.fromiter(
        (int(r[0]) for r in cn.execute("select primaryid from demo_latest order by rowid")),
        dtype=np.int64,
        count=n,
    )
    pids.sort()
    print("Sorted report index loaded; loading validated receipt dates", flush=True)
    ix = {}

    def idxget(pid):
        if pid in ix:
            return ix[pid]
        pos = int(np.searchsorted(pids, pid))
        answer = pos if pos < n and pids[pos] == pid else None
        ix[pid] = answer
        return answer

    arr = {
        "pids": pids,
        "ps": np.zeros(n, np.uint8),
        "psss": np.zeros(n, np.uint8),
        "pts": np.zeros(n, np.uint16),
        "fda": np.zeros(n, np.int32),
        "sus_first": np.zeros(n, np.uint64),
        "sus_multi": np.zeros(n, np.bool_),
        "sus_missing": np.zeros(n, np.bool_),
        "indi": np.zeros(n, np.uint8),
        "diarrhoea": np.zeros(n, np.bool_),
    }
    cursor = cn.execute("select primaryid,fda_dt from demo_latest")
    while batch := cursor.fetchmany(100000):
        indices = np.searchsorted(pids, np.array([int(r[0]) for r in batch], dtype=np.int64))
        arr["fda"][indices] = np.array([validdate(r[1]) for r in batch], dtype=np.int32)
    print("Building stringent DEMO duplicate-candidate blocks", flush=True)
    blocksql = """WITH eligible AS (SELECT primaryid,caseid,event_dt,init_fda_dt,age,age_cod,sex,UPPER(COALESCE(NULLIF(TRIM(occr_country),''),TRIM(reporter_country))) country FROM demo_latest WHERE length(TRIM(event_dt))=8 AND length(TRIM(init_fda_dt))=8 AND TRIM(age)<>'' AND age_cod IN ('YR','MON','WK','DY','DEC') AND sex IN ('F','M') AND COALESCE(NULLIF(TRIM(occr_country),''),TRIM(reporter_country))<>''), blocks AS (SELECT *,COUNT(*) OVER (PARTITION BY event_dt,init_fda_dt,CAST(age AS REAL),age_cod,sex,country) block_n FROM eligible) SELECT primaryid,caseid,event_dt,init_fda_dt,age,age_cod,sex,country FROM blocks WHERE block_n>1"""
    # SQLite groups on disk; retain only candidate-block records, avoiding a whole-database Python key dictionary.
    dupdemo = {}
    for pid, case, ev, init, age, unit, sex, country in cn.execute(blocksql):
        if not validdate(ev) or not validdate(init) or ageyears(age, unit) is None:
            continue
        key = (ev, init, str(float(age)), unit, sex, country)
        i = int(np.searchsorted(pids, int(pid)))
        dupdemo[i] = (int(case), key)
    print(f"Candidate DEMO block reports: {len(dupdemo):,}", flush=True)
    ddrug = defaultdict(set)
    dpt = defaultdict(set)
    dincomplete = set()
    target_seq = defaultdict(lambda: defaultdict(set))
    tindi = defaultdict(lambda: defaultdict(set))
    target_any = defaultdict(set)
    target_linkbits = defaultdict(lambda: defaultdict(int))
    tcats = defaultdict(lambda: [0, 0])
    tallnames = defaultdict(set)
    counts = Counter()
    sources = []
    zips = sorted(s5.FAERS_ZIPS.glob("faers_ascii_*.zip"))
    zips = [p for p in zips if "2019q1" <= p.stem[-6:].lower() <= "2026q1"]
    if len(zips) != 29:
        raise RuntimeError("Wrong source quarters")
    for qi, zp in enumerate(zips, 1):
        start = time.time()
        susseq = defaultdict(set)
        with zipfile.ZipFile(zp) as z:
            # DRUG first supplies exact suspect-sequence linkage for this release.
            for pid, c, h in rows(z, "DRUG", counts):
                i = idxget(pid)
                if i is None:
                    continue
                counts["DRUG_latest_rows"] += 1
                role = val(c, h, "role_cod").upper()
                ai = val(c, h, "prod_ai")
                name = val(c, h, "drugname")
                seq = val(c, h, "drug_seq")
                is_sus = role in ("PS", "SS")
                if is_sus:
                    if seq:
                        susseq[i].add(seq)
                    comps, ok = ingredients(ai)
                    if not ok:
                        arr["sus_missing"][i] = True
                    for comp in comps:
                        hcomp = hash64(comp)
                        if not arr["sus_first"][i]:
                            arr["sus_first"][i] = hcomp
                        elif arr["sus_first"][i] != hcomp:
                            arr["sus_multi"][i] = True
                if i in dupdemo:
                    comps, ok = ingredients(ai)
                    if not ok:
                        dincomplete.add(i)
                    else:
                        ddrug[i].update((role, x) for x in comps)
                # Match target using the same dictionary as the primary analysis.
                d = s5.match_drug(name, ai) if is_sus else None
                if d:
                    arr["psss"][i] |= 1 << DRUGS.index(d)
                if d and role == "PS":
                    arr["ps"][i] |= 1 << DRUGS.index(d)
                    if seq:
                        target_seq[i][d].add(seq)
                # Small target set may become known later in this quarter; store class bits for all records only if target known.
                # A second pass of only target DRUG rows is done below to make order independent.
            for pid, c, h in rows(z, "INDI", counts):
                i = idxget(pid)
                if i is None:
                    continue
                seq = val(c, h, "indi_drug_seq")
                pt = val(c, h, "indi_pt")
                bits = (
                    (1 if BREAST.search(pt) else 0)
                    | (2 if RENAL.search(pt) else 0)
                    | (4 if ONCO.search(pt) else 0)
                )
                if seq and seq in susseq.get(i, ()):
                    arr["indi"][i] |= bits
                if arr["ps"][i]:
                    for d, seqs in target_seq.get(i, {}).items():
                        if seq in seqs and pt:
                            tindi[i][d].add(pt)
                            target_linkbits[i][d] |= bits
            for pid, c, h in rows(z, "REAC", counts):
                i = idxget(pid)
                if i is None:
                    continue
                pt = val(c, h, "pt").strip().lower()
                if pt in PTBITS:
                    arr["pts"][i] |= PTBITS[pt]
                if pt in ("diarrhoea", "diarrhea"):
                    arr["diarrhoea"][i] = True
                if i in dupdemo and pt:
                    dpt[i].add(pt)
            # Same ZIP decompression repeated for DRUG only; avoids storing all concomitant text and report IDs in memory.
            for pid, c, h in rows(z, "DRUG", counts):
                i = idxget(pid)
                if i is None or not arr["ps"][i]:
                    continue
                ai = val(c, h, "prod_ai")
                name = val(c, h, "drugname")
                role = val(c, h, "role_cod").upper()
                tags = classtags(ai)
                tcats[i][0] |= tags
                if role in ("PS", "SS"):
                    tcats[i][1] |= tags
                tallnames[i].add((role, ai or name, bool(ai)))
        sources.append(
            {"file": str(zp), "bytes": zp.stat().st_size, "mtime_ns": zp.stat().st_mtime_ns}
        )
        print(
            f"{qi}/29 {zp.stem} done {time.time()-start:.1f}s; PS counts "
            + str({d: int(np.count_nonzero(arr["ps"] & (1 << k))) for k, d in enumerate(DRUGS)}),
            flush=True,
        )
        del susseq
        ix.clear()
    print("Completing duplicate candidates and target demographics", flush=True)
    dupgroups = defaultdict(list)
    for i, (case, key) in dupdemo.items():
        if i not in dincomplete and ddrug.get(i) and dpt.get(i):
            dupgroups[(key, tuple(sorted(ddrug[i])), tuple(sorted(dpt[i])))].append((case, i))
    remove = []
    dupsum = []
    for vals in dupgroups.values():
        if len(vals) < 2:
            continue
        vals.sort()
        remove.extend(i for case, i in vals[1:])
        dupsum.append(
            {
                "candidate_group_size": len(vals),
                "target_reports_n": sum(bool(arr["ps"][i]) for _, i in vals),
                "primary_reports_n": sum(bool(arr["pts"][i] & 7) for _, i in vals),
                "target_primary_reports_n": sum(
                    bool(arr["ps"][i] and arr["pts"][i] & 7) for _, i in vals
                ),
            }
        )
    arr["dup_remove"] = np.zeros(n, np.bool_)
    arr["dup_remove"][remove] = True
    targetidx = np.flatnonzero(arr["ps"])
    targetpids = [int(pids[i]) for i in targetidx]
    demo = {}
    cn.row_factory = sqlite3.Row
    for off in range(0, len(targetpids), 800):
        chunk = targetpids[off : off + 800]
        for r in cn.execute(
            "select * from demo_latest where primaryid in (" + ",".join("?" for _ in chunk) + ")",
            list(map(str, chunk)),
        ):
            demo[int(np.searchsorted(pids, int(r["primaryid"])))] = dict(r)
    cn.close()
    obj = {
        "demo": demo,
        "target_indications": {i: {d: sorted(v) for d, v in by.items()} for i, by in tindi.items()},
        "target_linkbits": {i: dict(by) for i, by in target_linkbits.items()},
        "target_classbits": dict(tcats),
        "target_drug_names": {i: sorted(v) for i, v in tallnames.items()},
        "duplicate_groups_aggregate": dupsum,
        "counts": dict(counts),
        "sources": sources,
        "duplicate_DEMO_candidate_reports_n": len(dupdemo),
        "duplicate_missing_drug_information_n": len(dincomplete),
    }
    np.savez_compressed(OUT / "local_report_cache_DO_NOT_PUBLISH.npz", **arr)
    with (OUT / "local_target_details_DO_NOT_PUBLISH.pkl").open("wb") as f:
        pickle.dump(obj, f, protocol=5)
    jwrite("source_manifest.json", sources)
    print("Extraction cached. Baseline audit will precede analysis.", flush=True)


def calc(label, drug, expo, cases, U, family):
    e = expo & U
    y = cases & U
    N = int(U.sum())
    a = int(np.count_nonzero(e & y))
    b = int(e.sum()) - a
    c = int(y.sum()) - a
    d = N - a - b - c
    if min(a, b, c, d) < 0:
        raise ValueError("Negative four-cell count")
    est = (a + b) > 0 and (c + d) > 0 and (a + c) > 0 and (b + d) > 0
    r = {
        "analysis": label,
        "family": family,
        "database": "FAERS",
        "drug": drug,
        "a": a,
        "b": b,
        "c": c,
        "d": d,
        "cohort_n": N,
        "target_n": a + b,
        "event_n": a + c,
        "estimable": est,
    }
    if est:
        ci = odds_ratio([[a, b], [c, d]], kind="conditional").confidence_interval()
        p = float(fisher_exact([[a, b], [c, d]]).pvalue)
        r.update(
            ROR=(a * d) / (b * c) if b * c else math.inf,
            exact_CI_low=float(ci.low),
            exact_CI_high=float(ci.high),
            Fisher_P=p,
        )
    else:
        r.update(ROR="", exact_CI_low="", exact_CI_high="", Fisher_P="")
    return r


def adjust(rs):
    for fam in sorted({x["family"] for x in rs}):
        g = [r for r in rs if r["family"] == fam and r["estimable"]]
        if g:
            for r, q in zip(g, multipletests([r["Fisher_P"] for r in g], method="fdr_by")[1]):
                r["BY_q"] = float(q)
                r["family_estimable_tests"] = len(g)
    for r in rs:
        r.setdefault("BY_q", "")
        r.setdefault("family_estimable_tests", 0)


def analyze():
    data = np.load(OUT / "local_report_cache_DO_NOT_PUBLISH.npz")
    a = {k: data[k] for k in data.files}
    data.close()
    with (OUT / "local_target_details_DO_NOT_PUBLISH.pkl").open("rb") as f:
        o = pickle.load(f)
    n = len(a["ps"])
    allU = np.ones(n, dtype=bool)
    primary = (a["pts"] & 7) > 0
    hyper = (a["pts"] & sum(1 << i for i in (10, 11, 12))) > 0
    exposure = {d: (a["ps"] & (1 << k)) > 0 for k, d in enumerate(DRUGS)}
    expectedN = [7630, 236, 2020, 14633, 203]
    expectedA = [162, 24, 72, 25, 0]
    audit = []
    for k, d in enumerate(DRUGS):
        nr = int(exposure[d].sum())
        nc = int(np.count_nonzero(exposure[d] & primary))
        ok = nr == expectedN[k] and nc == expectedA[k]
        audit.append(
            {
                "drug": d,
                "target_n": nr,
                "expected_target_n": expectedN[k],
                "primary_n": nc,
                "expected_primary_n": expectedA[k],
                "pass": ok,
            }
        )
    csvout("baseline_reconciliation.csv", audit)
    if not all(x["pass"] for x in audit):
        raise RuntimeError("STOP: baseline reconciliation failed; no further analysis")
    if n != 10463182 or int(primary.sum()) != 29665:
        raise RuntimeError(
            "STOP: full eligible report or primary-outcome background count differs from original analysis"
        )
    auditmeta = {
        "reports_n": n,
        "unique_pid_n": len(np.unique(a["pids"])),
        "primary_reports_all_drugs_n": int(primary.sum()),
        "missing_valid_FDA_DT_n": int(np.count_nonzero(a["fda"] == 0)),
        "sole_suspect_eligible_n": int(
            np.count_nonzero((a["sus_first"] > 0) & ~a["sus_multi"] & ~a["sus_missing"])
        ),
        "multiple_suspect_components_n": int(a["sus_multi"].sum()),
        "missing_or_ambiguous_suspect_ingredient_n": int(a["sus_missing"].sum()),
        "crossCASEID_candidate_groups_n": len(o["duplicate_groups_aggregate"]),
        "candidate_reports_removed_sensitivity_n": int(a["dup_remove"].sum()),
        "row_audit": o["counts"],
        "DRUG_processing_passes": 2,
        "DRUG_unique_source_lines": o["counts"].get("DRUG_raw_rows", 0) // 2,
        "DRUG_row_count_note": "DRUG_raw_rows includes TWO full processing passes; divide by two for source lines. INDI and REAC have one pass. Source lines are not unique reports.",
        "definitions_locked_before_extraction": True,
    }
    jwrite("audit_before_analysis.json", auditmeta)
    flow = [
        {
            "category": "no_PS_SS_ingredient_record",
            "reports_n": int(np.count_nonzero((a["sus_first"] == 0) & ~a["sus_missing"])),
        },
        {
            "category": "any_missing_or_ambiguous_PS_SS_ingredient",
            "reports_n": int(a["sus_missing"].sum()),
        },
        {
            "category": "multiple_distinct_components_without_missing",
            "reports_n": int(np.count_nonzero(a["sus_multi"] & ~a["sus_missing"])),
        },
        {
            "category": "one_unambiguous_encoded_component",
            "reports_n": int(
                np.count_nonzero((a["sus_first"] > 0) & ~a["sus_multi"] & ~a["sus_missing"])
            ),
        },
    ]
    if sum(r["reports_n"] for r in flow) != n:
        raise RuntimeError("Sole-suspect cohort audit categories do not reconcile")
    csvout("sole_suspect_cohort_flow.csv", flow)
    print(
        "Audit passes all five drug baselines; proceeding with locked exploratory definitions",
        flush=True,
    )
    desc = []
    indrows = []
    medrows = []
    miss = []
    component = []
    exclusive = []
    intersections = []
    sens = []
    temporal = []
    duplicates = []
    missingpatterns = []
    plausibility = []
    indstates = []
    drugterms = []
    ingredientchecks = []
    jwrite(
        "descriptive_missingness_dictionary.json",
        {
            "unknown_indication_terms": sorted(INDICATION_UNKNOWN_TERMS),
            "known_indication_rule": "At least one linked PT not in explicit unknown/unspecified indication dictionary. This does not validate the diagnosis.",
            "ingredient_completeness": "All listed drug records must have nonblank PROD_AI and pass the component parser; explicit unknown placeholders and ambiguous formulations are not known.",
            "co_record_note": "Computed from active-ingredient fields with curated dictionary; unknown or unmatched is not evidence of no use.",
        },
    )
    domains = {
        "primary": primary,
        "dehydration": (a["pts"] & 24) > 0,
        "electrolyte": (a["pts"] & 96) > 0,
        "renal": (a["pts"] & 896) > 0,
    }
    combo = np.zeros(n, np.uint8)
    for bit, mask in enumerate(domains.values()):
        combo |= mask.astype(np.uint8) << bit
    sole = (
        (a["sus_first"] > 0) & (a["sus_first"] != hash64("")) & ~a["sus_multi"] & ~a["sus_missing"]
    )
    normalized_empty_evidence = {
        "reports_first_suspect_normalizes_empty_n": int(
            np.count_nonzero(a["sus_first"] == hash64(""))
        ),
        "sole_candidate_normalizes_empty_n": int(
            np.count_nonzero((a["sus_first"] == hash64("")) & ~a["sus_multi"] & ~a["sus_missing"])
        ),
        "handling": "Such reports are excluded from sole-suspect restriction; raw main analysis remains unchanged.",
    }
    jwrite("normalized_empty_ingredient_audit.json", normalized_empty_evidence)
    for d in DRUGS:
        E = exposure[d]
        ids = np.flatnonzero(E & primary)
        ncase = len(ids)
        sex = Counter()
        age = Counter()
        ind = Counter()
        patterns = Counter()
        validages = []
        link_n = 0
        istates = Counter()
        iterms = Counter()
        for i0 in ids:
            i = int(i0)
            r = o["demo"][i]
            years = ageyears(r["age"], r["age_cod"])
            sx = str(r.get("sex") or "").strip().upper() or "blank"
            sex[sx] += 1
            group = (
                "unusable"
                if years is None
                else "<18"
                if years < 18
                else "18-44"
                if years < 45
                else "45-64"
                if years < 65
                else "65-74"
                if years < 75
                else ">=75"
            )
            age[group] += 1
            if years is not None:
                validages.append(years)
            pts = o["target_indications"].get(i, {}).get(d, [])
            isknown = any(known_indication(p) for p in pts)
            link_n += isknown
            istates["known"] += isknown
            istates["unknown_only"] += bool(pts) and not isknown
            istates["no_linked_record"] += not bool(pts)
            for p in pts:
                ind[p] += 1
            patterns[
                (
                    years is None,
                    sx not in ("F", "M"),
                    not isknown,
                    not bool(validdate(r["event_dt"])),
                )
            ] += 1
        for typ, cnt, levels in [
            ("age_group", age, ["<18", "18-44", "45-64", "65-74", ">=75", "unusable"]),
            ("sex_code", sex, sorted(set(sex) | {"F", "M", "blank"})),
        ]:
            for level in levels:
                desc.append(
                    {
                        "database": "FAERS",
                        "drug": d,
                        "primary_reports_n": ncase,
                        "variable": typ,
                        "category": level,
                        "n": cnt[level],
                        "pct_total_primary": 100 * cnt[level] / ncase if ncase else "",
                    }
                )
        desc.append(
            {
                "database": "FAERS",
                "drug": d,
                "primary_reports_n": ncase,
                "variable": "age_years",
                "category": "valid_median_Q1_Q3",
                "n": len(validages),
                "median": float(np.median(validages)) if validages else "",
                "Q1": float(np.quantile(validages, 0.25)) if validages else "",
                "Q3": float(np.quantile(validages, 0.75)) if validages else "",
            }
        )
        for p, cnt in ind.most_common():
            indrows.append(
                {
                    "drug": d,
                    "primary_reports_n": ncase,
                    "indication_PT": p,
                    "known_by_dictionary": known_indication(p),
                    "reports_n": cnt,
                    "pct_total_primary": 100 * cnt / ncase if ncase else "",
                    "nonexclusive": True,
                }
            )
        for state in ("known", "unknown_only", "no_linked_record"):
            indstates.append(
                {
                    "drug": d,
                    "primary_reports_n": ncase,
                    "indication_status": state,
                    "reports_n": istates[state],
                    "pct_primary": 100 * istates[state] / ncase if ncase else "",
                }
            )
        for var, available in [
            ("age_valid", len(validages)),
            ("sex_F_M", sex["F"] + sex["M"]),
            ("target_linked_indication", link_n),
            ("valid_event_date", sum(bool(validdate(o["demo"][int(i)]["event_dt"])) for i in ids)),
        ]:
            miss.append(
                {
                    "drug": d,
                    "variable": var,
                    "primary_reports_n": ncase,
                    "available_n": available,
                    "missing_or_unusable_n": ncase - available,
                    "missing_or_unusable_pct": 100 * (ncase - available) / ncase if ncase else "",
                }
            )
        missing_ai = sum(
            any(
                not row[2] or not ingredients(row[1])[1]
                for row in o["target_drug_names"].get(int(i), [])
            )
            for i in ids
        )
        miss.append(
            {
                "drug": d,
                "variable": "PROD_AI_known_unambiguous_for_all_listed_drug_rows",
                "primary_reports_n": ncase,
                "available_n": ncase - missing_ai,
                "missing_or_unusable_n": missing_ai,
                "missing_or_unusable_pct": 100 * missing_ai / ncase if ncase else "",
            }
        )
        future_events = 0
        negative_age = 0
        over120_age = 0
        for i in ids:
            row = o["demo"][int(i)]
            ev = validdate(row["event_dt"])
            fd = validdate(row["fda_dt"])
            future_events += bool(ev and fd and ev > fd)
            try:
                yr = float(row["age"]) * AGE_FACTORS.get(str(row["age_cod"]).upper(), math.nan)
                negative_age += yr < 0
                over120_age += yr > 120
            except (ValueError, TypeError):
                pass
        plausibility.append(
            {
                "drug": d,
                "primary_reports_n": ncase,
                "event_date_after_latest_FDA_date_n": future_events,
                "negative_age_n": negative_age,
                "age_above_120_n": over120_age,
                "action": "flagged; no automatic deletion or correction of source reports",
            }
        )
        for pat, cnt in patterns.items():
            missingpatterns.append(
                {
                    "drug": d,
                    "age_unusable": pat[0],
                    "sex_unknown": pat[1],
                    "indication_not_known": pat[2],
                    "event_date_unusable": pat[3],
                    "reports_n": cnt,
                }
            )
        caseclass = {}
        for i0 in ids:
            i = int(i0)
            bits = [0, 0]
            for role, ai, present in o["target_drug_names"].get(i, []):
                if not present:
                    continue
                tag = classtags(ai)
                bits[0] |= tag
                if role in ("PS", "SS"):
                    bits[1] |= tag
                if tag:
                    iterms[(role, ai, tag)] += 1
            caseclass[i] = bits
        for bit, cat in enumerate(CATS):
            for which, rolelabel in [(0, "any_reported_role"), (1, "PS_or_SS")]:
                count = sum(bool(caseclass[int(i)][which] & (1 << bit)) for i in ids)
                medrows.append(
                    {
                        "drug": d,
                        "primary_reports_n": ncase,
                        "ingredient_class": cat,
                        "role_scope": rolelabel,
                        "co_record_n": count,
                        "pct_primary": 100 * count / ncase if ncase else "",
                    }
                )
        for (role, term, tag), count in iterms.most_common():
            drugterms.append(
                {
                    "drug": d,
                    "reported_role": role,
                    "reported_PROD_AI": term,
                    "dictionary_classes": ";".join(
                        name for bit, name in enumerate(CATS) if tag & (1 << bit)
                    ),
                    "primary_reports_n_with_term": count,
                }
            )
        broad = E & (combo > 0)
        nb = int(broad.sum())
        for name, mask in domains.items():
            ct = int(np.count_nonzero(E & mask))
            component.append(
                {
                    "database": "FAERS",
                    "drug": d,
                    "domain": name,
                    "domain_n": ct,
                    "broad_union_n": nb,
                    "pct_broad_nonexclusive": 100 * ct / nb if nb else "",
                }
            )
        for code in range(1, 16):
            exclusive.append(
                {
                    "database": "FAERS",
                    "drug": d,
                    "combination_code": code,
                    "exclusive_domains": "+".join(
                        name for bit, name in enumerate(domains) if code & (1 << bit)
                    ),
                    "n": int(np.count_nonzero(E & (combo == code))),
                    "broad_union_n": nb,
                }
            )
        for size in (2, 3, 4):
            for names in itertools.combinations(domains, size):
                mask = E.copy()
                for name in names:
                    mask &= domains[name]
                intersections.append(
                    {
                        "database": "FAERS",
                        "drug": d,
                        "intersection": " & ".join(names),
                        "n": int(mask.sum()),
                        "broad_union_n": nb,
                    }
                )
        duplicates.append(
            {
                "drug": d,
                "target_reports_n": int(E.sum()),
                "target_candidate_removed_n": int(np.count_nonzero(E & a["dup_remove"])),
                "primary_reports_n": ncase,
                "primary_candidate_removed_n": int(np.count_nonzero(E & primary & a["dup_remove"])),
            }
        )
        inconsistent = E & sole & (a["sus_first"] != hash64(d.upper()))
        ingredientchecks.append(
            {
                "drug": d,
                "sole_suspect_target_candidate_n": int(np.count_nonzero(E & sole)),
                "sole_suspect_target_primary_candidate_n": int(
                    np.count_nonzero(E & sole & primary)
                ),
                "sole_component_not_matching_target_n": int(inconsistent.sum()),
                "sole_component_not_matching_target_primary_n": int(
                    np.count_nonzero(inconsistent & primary)
                ),
                "action": "Remove inconsistent target reports from this sensitivity universe; do not reclassify them as comparator reports.",
            }
        )
        sens.append(
            calc(
                "sole_known_suspect_active_ingredient",
                d,
                E,
                primary,
                sole & ~inconsistent,
                "stage8_non_temporal_sensitivity",
            )
        )
        sens.append(
            calc(
                "candidate_duplicate_removal",
                d,
                E,
                primary,
                ~a["dup_remove"],
                "stage8_non_temporal_sensitivity",
            )
        )
        for cat, bit in [("breast", 1), ("renal_cancer", 2)]:
            targetlinked = np.zeros(n, bool)
            for i, by in o["target_linkbits"].items():
                if by.get(d, 0) & bit:
                    targetlinked[i] = True
            U = ((a["indi"] & bit) > 0) & ~(E & ~targetlinked)
            sens.append(
                calc(
                    "exact_suspect_linked_" + cat,
                    d,
                    targetlinked,
                    primary,
                    U,
                    "stage8_non_temporal_sensitivity",
                )
            )
        for name, mask in [("hyperglycemia", hyper), ("diarrhoea", a["diarrhoea"])]:
            sens.append(
                calc(
                    "symmetric_exclusion_all_" + name + "_coreports",
                    d,
                    E,
                    primary,
                    ~mask,
                    "stage8_non_temporal_sensitivity",
                )
            )
        if d in ("everolimus", "temsirolimus"):
            competing = (a["psss"] & 7) > 0
            sens.append(
                calc(
                    "symmetric_exclusion_PI3K_AKT_PSSS_coreports",
                    d,
                    E,
                    primary,
                    ~competing,
                    "stage8_non_temporal_sensitivity",
                )
            )
        for name, lo, hi in [
            ("common_complete_2025", 20250101, 20251231),
            ("common_complete_2024_2025", 20240101, 20251231),
        ]:
            U = (a["fda"] >= lo) & (a["fda"] <= hi)
            sens.append(calc(name, d, E, primary, U, "stage8_non_temporal_sensitivity"))
    years = a["fda"] // 10000
    months = (a["fda"] // 100) % 100
    quarters = years * 10 + (months - 1) // 3 + 1
    for unit, strata in [("year", years), ("quarter", quarters)]:
        for code in np.unique(strata[a["fda"] > 0]).tolist():
            U = (strata == code) & (a["fda"] > 0)
            for d in DRUGS:
                r = calc(
                    f"FDA_DT_{unit}_{code}",
                    d,
                    exposure[d],
                    primary,
                    U,
                    "stage8_all_" + unit + "_strata",
                )
                r.update(
                    time_unit=unit,
                    time_value=code,
                    source_coverage="partial_2026Q1"
                    if code == 2026 or code == 20261
                    else "pre2019_receipt_in_later_release"
                    if (code if unit == "year" else code // 10) < 2019
                    else "within_source_calendar",
                )
                temporal.append(r)
    adjust(sens)
    adjust(temporal)
    for r in sens:
        k = DRUGS.index(r["drug"])
        r.update(
            reference_cohort_n=n,
            reference_target_n=expectedN[k],
            reference_target_primary_n=expectedA[k],
            excluded_total_reports_n=n - r["cohort_n"],
            excluded_target_reports_n=expectedN[k] - r["target_n"],
            excluded_target_primary_reports_n=expectedA[k] - r["a"],
            excluded_comparator_reports_n=(n - expectedN[k]) - (r["c"] + r["d"]),
        )
    for name, rs in [
        ("primary_case_demographic_distribution.csv", desc),
        ("primary_case_linked_indications.csv", indrows),
        ("primary_case_indication_status.csv", indstates),
        ("primary_case_missingness.csv", miss),
        ("primary_case_missingness_patterns.csv", missingpatterns),
        ("primary_case_plausibility_flags.csv", plausibility),
        ("primary_case_drug_class_corecords.csv", medrows),
        ("primary_case_drug_dictionary_matches.csv", drugterms),
        ("sole_suspect_target_ingredient_consistency.csv", ingredientchecks),
        ("broad_domain_contribution.csv", component),
        ("broad_15_exclusive_combinations.csv", exclusive),
        ("broad_all_domain_intersections.csv", intersections),
        ("sensitivity_exact_ROR_BY.csv", sens),
        ("all_year_quarter_exact_ROR_BY.csv", temporal),
        ("crossCASEID_candidate_impact.csv", duplicates),
        ("crossCASEID_candidate_group_aggregates.csv", o["duplicate_groups_aggregate"]),
    ]:
        csvout(name, rs)
    validation = []
    for d in DRUGS:
        total = sum(r["n"] for r in exclusive if r["drug"] == d)
        broad = next(r["broad_union_n"] for r in component if r["drug"] == d)
        validation.append(
            {
                "check": "15 mutually exclusive cells sum to broad union",
                "drug": d,
                "observed": total,
                "expected": broad,
                "pass": total == broad,
            }
        )
        for typ in ("age_group", "sex_code"):
            vals = [r for r in desc if r["drug"] == d and r["variable"] == typ]
            observed = sum(r["n"] for r in vals)
            expected = expectedA[DRUGS.index(d)]
            validation.append(
                {
                    "check": typ + " including missing sums to cases",
                    "drug": d,
                    "observed": observed,
                    "expected": expected,
                    "pass": observed == expected,
                }
            )
        status_sum = sum(r["reports_n"] for r in indstates if r["drug"] == d)
        validation.append(
            {
                "check": "known+unknown_only+no_linked_indication sums to cases",
                "drug": d,
                "observed": status_sum,
                "expected": expectedA[DRUGS.index(d)],
                "pass": status_sum == expectedA[DRUGS.index(d)],
            }
        )
    csvout("analysis_validation.csv", validation)
    if not all(r["pass"] for r in validation):
        raise RuntimeError("Final reconciliation failed")
    jwrite(
        "analysis_completion.json",
        {
            "completed_at": datetime.now().isoformat(timespec="seconds"),
            "audit_pass": True,
            "validation_pass": True,
            "new_sensitivity_tests": len(sens),
            "new_sensitivity_estimable": sum(r["estimable"] for r in sens),
            "temporal_tests": len(temporal),
            "temporal_estimable": sum(r["estimable"] for r in temporal),
            "cautions": [
                "All new analyses are post-review exploratory.",
                "Candidate duplicates are not proven duplicate clinical cases.",
                "Symmetric co-report/background restriction is not a validated masking correction and may introduce selection bias.",
                "Public outputs are aggregate only; local NPZ/PKL caches must not be uploaded.",
            ],
        },
    )
    print("COMPLETE", flush=True)


if __name__ == "__main__":
    if "--extract" in sys.argv:
        extract()
    elif "--analyze" in sys.argv:
        analyze()
    else:
        if not (OUT / "local_report_cache_DO_NOT_PUBLISH.npz").exists():
            extract()
        analyze()
