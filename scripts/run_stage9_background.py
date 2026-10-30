"""Revision-stage FAERS drug-scope and competing-background audit.

Audit and analysis are separate commands. Report-level caches and candidate rows
are local only. No clinical sequence, causal effect or incidence is inferred.
"""
from __future__ import annotations
import argparse, csv, gzip, hashlib, importlib.util, json, math, re, time, zipfile
from collections import Counter
from datetime import datetime
from pathlib import Path
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
from statsmodels.stats.multitest import multipletests

from package_paths import OUTPUT, PACKAGE

OUT = OUTPUT / "03_确认后补充分析/阶段九_选药范围与竞争背景"
ROOT = OUTPUT
S8 = OUTPUT / "03_确认后补充分析/阶段八_核验后修复/FAERS"
S5 = Path(__file__).resolve().parent / "run_stage5_primary_analysis.py"
spec = importlib.util.spec_from_file_location("s5", S5)
s5 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s5)
DRUGS = list(s5.DRUG_TERMS)
SG = {
    "canagliflozin": ["CANAGLIFLOZIN", "INVOKANA", "INVOKAMET", "VOKANAMET", "CANAGLU"],
    "dapagliflozin": ["DAPAGLIFLOZIN", "FARXIGA", "FORXIGA", "XIGDUO", "QTERN", "QTERNMET"],
    "empagliflozin": ["EMPAGLIFLOZIN", "JARDIANCE", "GLYXAMBI", "SYNJARDY", "TRIJARDY"],
    "ertugliflozin": ["ERTUGLIFLOZIN", "STEGLATRO", "STEGLUJAN", "SEGLUROMET"],
    "bexagliflozin": ["BEXAGLIFLOZIN", "BRENZAVVY"],
    "sotagliflozin": ["SOTAGLIFLOZIN", "INPEFA", "ZYNQUISTA"],
    "ipragliflozin": ["IPRAGLIFLOZIN", "SUGLAT"],
    "luseogliflozin": ["LUSEOGLIFLOZIN", "LUSEFI"],
    "tofogliflozin": ["TOFOGLIFLOZIN", "DEBERZA", "APLEWAY"],
}
SCOPE = {
    "nab_sirolimus_explicit": [
        "FYARRO",
        "NAB SIROLIMUS",
        "SIROLIMUS PROTEIN BOUND",
        "ALBUMIN BOUND SIROLIMUS",
    ],
    "sirolimus_any_formulation": ["SIROLIMUS", "RAPAMUNE", "FYARRO"],
    "idelalisib": ["IDELALISIB", "ZYDELIG"],
    "duvelisib": ["DUVELISIB", "COPIKTRA"],
    "copanlisib": ["COPANLISIB", "ALIQOPA"],
    "umbralisib": ["UMBRALISIB", "UKONIQ"],
    "leniolisib": ["LENIOLISIB", "JOENJA"],
}
ALIAS_SOURCES = [
    "https://www.fda.gov/drugs/fda-adverse-event-monitoring-system-aems/july-september-2023-potential-signals-serious-risksnew-safety-information-identified-fda-adverse",
    "https://www.pmda.go.jp/files/000208624.pdf",
    "https://www.pmda.go.jp/files/000197771.pdf",
    "https://www.ema.europa.eu/en/medicines/human/EPAR/zynquista",
    "https://www.ema.europa.eu/en/medicines/human/EPAR/vokanamet",
    "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2022/213312Orig1s000Approv.pdf",
]


def jwrite(name, obj):
    (OUT / name).write_text(
        json.dumps(obj, ensure_ascii=False, indent=2, default=str), encoding="utf-8"
    )


def readcsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def writecsv(name, rows):
    if not rows:
        return
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for row in rows for k in row)))
        w.writeheader()
        w.writerows(rows)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def norm(s):
    return re.sub("[^A-Z0-9]+", " ", s.upper()).strip()


def compile_terms(mapping):
    return {
        k: re.compile(
            r"(?<![A-Z0-9])(?:"
            + "|".join(re.escape(norm(v)).replace(r"\ ", r"\s*") for v in vals)
            + r")(?![A-Z])"
        )
        for k, vals in mapping.items()
    }


INITIAL_SG = {k: list(v) for k, v in SG.items()}
FINAL_DICTIONARY = PACKAGE / "definitions" / "sglt2_dictionary.json"
SG = json.loads(FINAL_DICTIONARY.read_text(encoding="utf-8"))["final_dictionary"]
SG_PAT = compile_terms(SG)
SCOPE_PAT = compile_terms(SCOPE)
QUICK = re.compile(
    "GLIFLOZIN|SIROLIMUS|FYARRO|IDELALISIB|DUVELISIB|COPANLISIB|UMBRALISIB|LENIOLISIB|"
    + "|".join(
        re.escape(x)
        for x in sorted(
            set(
                v for vals in list(SG.values()) + list(SCOPE.values()) for v in vals if " " not in v
            )
        )
    )
)
BLANK = {
    "",
    "UNKNOWN",
    "UNK",
    "UNSPECIFIED",
    "NOT REPORTED",
    "NOT PROVIDED",
    "NOT KNOWN",
    "NA",
    "N A",
    "NULL",
    "NONE",
}


def lock_plan():
    p = PACKAGE / "definitions" / "background_analysis_rules.json"
    return json.loads(p.read_text(encoding="utf-8"))


def audit():
    lock_plan()
    t0 = time.time()
    cache = S8 / "local_report_cache_DO_NOT_PUBLISH.npz"
    with np.load(cache) as z:
        pids = z["pids"]
    assert len(pids) == 10463182 and np.all(pids[1:] > pids[:-1])
    n = len(pids)
    arrays = {
        x: np.zeros(n, np.uint16)
        for x in ["sg_ps", "sg_psss", "sg_any", "scope_ps", "scope_psss", "scope_any"]
    }
    arrays["field_flags"] = np.zeros(n, np.uint8)
    # field flags: present DRUG, any missing name, any missing ingredient,
    # any unknown/missing role, any missing drug sequence.
    counts = Counter()
    terms = Counter()
    unknown = Counter()
    sources = []
    duplicate_candidates = Counter()
    alias_match_cache = {}
    candidate_seen = set()
    candidate_path = OUT / "local_candidate_drug_rows_DO_NOT_PUBLISH.csv.gz"
    zips = sorted(s5.FAERS_ZIPS.glob("faers_ascii_*.zip"))
    zips = [p for p in zips if "2019q1" <= p.stem[-6:].lower() <= "2026q1"]
    assert len(zips) == 29
    previous = {
        x["filename"]: {"bytes": int(x["bytes"]), "sha256": x["sha256"]}
        for x in readcsv(PACKAGE / "definitions/source_files.csv") if x["database"] == "FAERS"
    }

    def match(name, ai):
        key = (name, ai)
        if key in alias_match_cache:
            return alias_match_cache[key]
        fields = (norm(name), norm(ai))
        sg = 0
        sc = 0
        for k, pat in enumerate(SG_PAT.values()):
            if any(pat.search(v) for v in fields):
                sg |= 1 << k
        for k, pat in enumerate(SCOPE_PAT.values()):
            if any(pat.search(v) for v in fields):
                sc |= 1 << k
        unknown_glif = "GLIFLOZIN" in " ".join(fields) and not sg
        ans = (sg, sc, unknown_glif)
        alias_match_cache[key] = ans
        return ans

    with gzip.open(candidate_path, "wt", encoding="utf-8", newline="") as cf:
        cw = csv.writer(cf)
        cw.writerow(
            [
                "PRIMARYID",
                "DRUG_SEQ",
                "role",
                "drugname",
                "PROD_AI",
                "sg_bits",
                "scope_bits",
                "unmapped_gliflozin",
            ]
        )
        for qi, zp in enumerate(zips, 1):
            start = time.time()
            stat = zp.stat()
            old = previous[zp.name]
            assert stat.st_size == old["bytes"], f"Source size changed {zp}"
            # Hash source ZIPs, not merely the filename; one recorded digest per release.
            digest = sha(zp)
            assert digest == old["sha256"], f"Source checksum changed {zp}"
            with zipfile.ZipFile(zp) as zf:
                member = s5.find_member(zf, "DRUG")
                assert member
                with zf.open(member) as f:
                    header = f.readline().decode("latin-1").rstrip("\r\n").lower().split("$")
                    h = {v: k for k, v in enumerate(header)}
                    assert all(
                        x in h for x in ["primaryid", "drug_seq", "role_cod", "drugname", "prod_ai"]
                    )
                    batch = []
                    flags = []
                    candidates = []

                    def flush():
                        if not batch:
                            return
                        rawids = np.asarray(batch, dtype=np.int64)
                        idx = np.searchsorted(pids, rawids)
                        valid = idx < n
                        valid[valid] &= pids[idx[valid]] == rawids[valid]
                        positions = idx[valid]
                        np.bitwise_or.at(
                            arrays["field_flags"],
                            positions,
                            np.asarray(flags, dtype=np.uint8)[valid],
                        )
                        counts["DRUG_latest_rows"] += int(valid.sum())
                        for loc, seq, role, name, ai, sg, sc, unmapped in candidates:
                            if not valid[loc]:
                                continue
                            i = int(idx[loc])
                            pid = int(rawids[loc])
                            key = (pid, seq, role, name, ai)
                            if key in candidate_seen:
                                counts["exact_duplicate_matched_rows"] += 1
                            else:
                                candidate_seen.add(key)
                            arrays["sg_any"][i] |= sg
                            arrays["scope_any"][i] |= sc
                            if role in ("PS", "SS"):
                                arrays["sg_psss"][i] |= sg
                                arrays["scope_psss"][i] |= sc
                            if role == "PS":
                                arrays["sg_ps"][i] |= sg
                                arrays["scope_ps"][i] |= sc
                            terms[(name, ai, role, sg, sc, unmapped)] += 1
                            if unmapped:
                                unknown[(name, ai, role)] += 1
                            cw.writerow([pid, seq, role, name, ai, sg, sc, int(unmapped)])
                        batch.clear()
                        flags.clear()
                        candidates.clear()

                    for ln in f:
                        counts["DRUG_source_rows"] += 1
                        c = ln.rstrip(b"\r\n").split(b"$")
                        if len(c) < len(header):
                            counts["short_rows"] += 1
                            continue

                        def v(k):
                            return c[h[k]].strip(b'" ').decode("latin-1")

                        try:
                            pid = int(v("primaryid"))
                        except ValueError:
                            counts["invalid_primaryid"] += 1
                            continue
                        name = v("drugname")
                        ai = v("prod_ai")
                        role = v("role_cod").upper()
                        seq = v("drug_seq")
                        fl = (
                            1
                            | (2 if name.upper() in BLANK else 0)
                            | (4 if ai.upper() in BLANK else 0)
                            | (8 if role not in ("PS", "SS", "C", "I") else 0)
                            | (16 if not seq else 0)
                        )
                        loc = len(batch)
                        batch.append(pid)
                        flags.append(fl)
                        if QUICK.search(name.upper()) or QUICK.search(ai.upper()):
                            sg, sc, unmapped = match(name, ai)
                            if sg or sc or unmapped:
                                candidates.append((loc, seq, role, name, ai, sg, sc, unmapped))
                        if len(batch) >= 100000:
                            flush()
                    flush()
            sources.append(
                {
                    "release": zp.stem,
                    "bytes": stat.st_size,
                    "sha256": digest,
                    "drug_member": member,
                    "drug_columns_n": len(header),
                    "drug_columns": header,
                }
            )
            print(
                f'{qi}/29 {zp.stem}: {time.time()-start:.1f}s; latest rows {counts["DRUG_latest_rows"]:,}',
                flush=True,
            )
    np.savez_compressed(OUT / "local_background_cache_DO_NOT_PUBLISH.npz", **arrays)
    writecsv(
        "matched_drug_terms_audit.csv",
        [
            dict(
                drugname=k[0],
                PROD_AI=k[1],
                role=k[2],
                sg_bits=k[3],
                scope_bits=k[4],
                unmapped_gliflozin=k[5],
                source_row_n=v,
            )
            for k, v in sorted(terms.items())
        ],
    )
    writecsv(
        "unmapped_gliflozin_terms.csv",
        [
            dict(drugname=k[0], PROD_AI=k[1], role=k[2], source_row_n=v)
            for k, v in sorted(unknown.items())
        ]
        or [{"drugname": "NONE", "PROD_AI": "", "role": "", "source_row_n": 0}],
    )
    field = arrays["field_flags"]
    missing = []
    for bit, var in [
        (1, "DRUG_record"),
        (2, "drugname"),
        (4, "PROD_AI"),
        (8, "role_cod"),
        (16, "DRUG_SEQ"),
    ]:
        num = (
            int(np.count_nonzero((field & 1) == 0))
            if bit == 1
            else int(np.count_nonzero(field & bit))
        )
        missing.append(
            {
                "variable": var,
                "reports_n": n,
                "missing_or_any_missing_record_n": num,
                "pct_reports": 100 * num / n,
                "definition": "No drug row"
                if bit == 1
                else "At least one DRUG row with missing/unknown field; other rows can contain usable information",
            }
        )
    writecsv("drug_fields_missingness.csv", missing)
    vals, cnt = np.unique(field, return_counts=True)
    writecsv(
        "drug_fields_missingness_patterns.csv",
        [
            {
                "flags": int(v),
                "reports_n": int(c),
                "has_drug": bool(v & 1),
                "any_missing_name": bool(v & 2),
                "any_missing_ingredient": bool(v & 4),
                "any_unknown_role": bool(v & 8),
                "any_missing_sequence": bool(v & 16),
            }
            for v, c in zip(vals, cnt)
        ],
    )
    jwrite("source_manifest.json", sources)
    jwrite(
        "audit_extraction.json",
        {
            "completed_at": datetime.now().isoformat(timespec="seconds"),
            "runtime_seconds": time.time() - t0,
            "reports_n": n,
            "unique_sorted_PRIMARYID": True,
            "counts": dict(counts),
            "sglt2_candidate_reports_any_role_n": int(np.count_nonzero(arrays["sg_any"])),
            "unmapped_name_role_combinations_n": len(unknown),
            "index_cache_sha256": sha(cache),
            "candidate_rows_private": True,
            "new_outcomes_or_effects_not_yet_summarized": True,
            "analysis_allowed": False,
        },
    )
    print(
        "Extraction finished. Audit review and an explicit audit acceptance file are required before analysis.",
        flush=True,
    )


def calc(label, drug, expo, Y, U):
    e = expo & U
    y = Y & U
    N = int(U.sum())
    a = int(np.count_nonzero(e & y))
    b = int(e.sum()) - a
    c = int(y.sum()) - a
    d = N - a - b - c
    ok = (a + b) > 0 and (c + d) > 0 and (a + c) > 0 and (b + d) > 0
    r = {
        "analysis": label,
        "drug": drug,
        "a": a,
        "b": b,
        "c": c,
        "d": d,
        "eligible_reports_n": N,
        "target_reports_n": a + b,
        "primary_reports_n": a + c,
        "estimable": ok,
    }
    if ok:
        ci = odds_ratio([[a, b], [c, d]], kind="conditional").confidence_interval()
        R = (a * d) / (b * c) if b * c else math.inf
        r.update(
            ROR=R,
            exact_CI_low=float(ci.low),
            exact_CI_high=float(ci.high),
            Fisher_P=float(fisher_exact([[a, b], [c, d]]).pvalue),
        )
    else:
        r.update(ROR="", exact_CI_low="", exact_CI_high="", Fisher_P="")
    return r


def analyze():
    plan = lock_plan()
    gate = json.loads((OUT / "audit_acceptance.json").read_text(encoding="utf-8"))
    assert gate["analysis_allowed"] and gate["unresolved_candidate_decision_n"] == 0
    assert gate["dictionary_amendment_sha256"] == sha(FINAL_DICTIONARY)
    with np.load(S8 / "local_report_cache_DO_NOT_PUBLISH.npz") as z:
        base = {k: z[k] for k in ["pids", "ps", "psss", "pts", "fda"]}
    with np.load(OUT / "local_background_cache_DO_NOT_PUBLISH.npz") as z:
        bg = {k: z[k] for k in z.files}
    n = len(base["pids"])
    U = np.ones(n, bool)
    Y = (base["pts"] & 7) > 0
    DKA = (base["pts"] & 1) > 0
    assert n == 10463182 and int(Y.sum()) == 29665
    expectedN = [7630, 236, 2020, 14633, 203]
    expectedA = [162, 24, 72, 25, 0]
    baselines = []
    for k, d in enumerate(DRUGS):
        E = (base["ps"] & (1 << k)) > 0
        r = calc("baseline", d, E, Y, U)
        assert r["a"] == expectedA[k] and r["target_reports_n"] == expectedN[k]
        baselines.append(r)
    writecsv("baseline_reconciliation.csv", baselines)
    # Data-audit and descriptive summaries precede new effect estimation.
    desc = []
    co = []
    scope = []
    for role in ["ps", "psss", "any"]:
        for k, name in enumerate(SG):
            m = (bg["sg_" + role] & (1 << k)) > 0
            desc.append(
                {
                    "ingredient": name,
                    "role": role,
                    "reports_n": int(m.sum()),
                    "primary_outcome_reports_n": int(np.count_nonzero(m & Y)),
                    "DKA_PT_reports_n": int(np.count_nonzero(m & DKA)),
                    "pct_all_primary_reports": 100 * np.count_nonzero(m & Y) / Y.sum(),
                }
            )
        m = bg["sg_" + role] > 0
        desc.append(
            {
                "ingredient": "ANY_INCLUDED_SGLT2_ACTIVE",
                "role": role,
                "reports_n": int(m.sum()),
                "primary_outcome_reports_n": int(np.count_nonzero(m & Y)),
                "DKA_PT_reports_n": int(np.count_nonzero(m & DKA)),
                "pct_all_primary_reports": 100 * np.count_nonzero(m & Y) / Y.sum(),
            }
        )
        for k, d in enumerate(DRUGS):
            E = (base["ps"] & (1 << k)) > 0
            counts = {
                "drug": d,
                "sglt2_role": role,
                "target_reports_n": int(E.sum()),
                "target_primary_reports_n": int(np.count_nonzero(E & Y)),
                "sglt2_coreports_n": int(np.count_nonzero(E & m)),
                "sglt2_primary_coreports_n": int(np.count_nonzero(E & m & Y)),
            }
            co.append(counts)
        for k, name in enumerate(SCOPE):
            m = (bg["scope_" + role] & (1 << k)) > 0
            scope.append(
                {
                    "candidate": name,
                    "role": role,
                    "reports_n": int(m.sum()),
                    "primary_outcome_reports_n": int(np.count_nonzero(m & Y)),
                    "DKA_PT_reports_n": int(np.count_nonzero(m & DKA)),
                }
            )
    writecsv("sglt2_background_contributions.csv", desc)
    writecsv("sglt2_target_coreports.csv", co)
    writecsv("candidate_scope_counts.csv", scope)
    results = []
    for k, d in enumerate(DRUGS):
        if d not in ("everolimus", "temsirolimus"):
            continue
        E = (base["ps"] & (1 << k)) > 0
        original = baselines[k]
        for role in ["psss", "any"]:
            exclusion = bg["sg_" + role] > 0
            r = calc("exclude_SGLT2_" + role, d, E, Y, ~exclusion)
            for cell in ["a", "b", "c", "d"]:
                r["original_" + cell] = original[cell]
                r["excluded_" + cell] = original[cell] - r[cell]
            r["original_ROR"] = original["ROR"]
            r["ROR_after_before_ratio"] = (
                r["ROR"] / original["ROR"]
                if r["estimable"]
                and original["ROR"] not in ("", 0)
                and math.isfinite(original["ROR"])
                else ""
            )
            r["excluded_background_primary_pct"] = 100 * r["excluded_c"] / original["c"]
            r[
                "interpretation"
            ] = "Descriptive reporting-background sensitivity; not a causal masking correction"
            assert all(r["excluded_" + c] >= 0 for c in ["a", "b", "c", "d"])
            results.append(r)
    for role in ["ps", "psss"]:
        E = (bg["scope_" + role] & 1) > 0
        r = calc("explicit_nab_sirolimus_" + role, "nab-sirolimus (explicit formulation)", E, Y, U)
        r[
            "interpretation"
        ] = "Supplementary formulation-specific FAERS estimate; not added to original primary family"
        results.append(r)
    valid = [r for r in results if r["estimable"]]
    qs = multipletests([r["Fisher_P"] for r in valid], method="fdr_by")[1] if valid else []
    for r, q in zip(valid, qs):
        r["BY_q"] = float(q)
    for r in results:
        r.setdefault("BY_q", "")
        r["family"] = "stage9_six_planned_contrasts"
        r["family_estimable_tests"] = len(valid)
    writecsv("additional_sensitivity_exact_BY.csv", results)
    jwrite(
        "analysis_completion.json",
        {
            "completed_at": datetime.now().isoformat(timespec="seconds"),
            "reports_n": n,
            "planned_contrasts": 6,
            "estimable_contrasts": len(valid),
            "all_baselines_match": True,
            "audit_gate": gate,
            "plan_sha256": sha(PACKAGE / "definitions" / "background_analysis_rules.json"),
            "dictionary_amendment_sha256": sha(FINAL_DICTIONARY),
            "no_changes_to_original_tests": True,
        },
    )
    for r in results:
        print(r, flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["plan", "audit", "analyze"])
    args = ap.parse_args()
    if args.action == "plan":
        lock_plan()
    elif args.action == "audit":
        audit()
    else:
        analyze()
