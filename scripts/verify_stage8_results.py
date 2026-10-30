"""Independent aggregate consistency and multiplicity checks; no source rescan."""
from pathlib import Path
import csv, json, math
import numpy as np

from package_paths import OUTPUT

OUT = OUTPUT / "03_确认后补充分析" / "阶段八_核验后修复" / "FAERS"
STAGE6 = OUTPUT / "03_确认后补充分析" / "阶段六_敏感性分析"


def read(p):
    with p.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write(name, rows):
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        w.writeheader()
        w.writerows(rows)


old = read(STAGE6 / "阶段六_敏感性ROR结果.csv")
oldyears = {
    (r["drug"], int(float(r["report_year_by_FDA_DT"]))): r
    for r in old
    if r.get("report_year_by_FDA_DT")
}
newtime = read(OUT / "all_year_quarter_exact_ROR_BY.csv")
newyears = {(r["drug"], int(r["time_value"])): r for r in newtime if r["time_unit"] == "year"}
rows = []
for key in sorted(set(oldyears) | set(newyears)):
    a = oldyears.get(key)
    b = newyears.get(key)
    r = {
        "drug": key[0],
        "year": key[1],
        "present_in_old": bool(a),
        "present_in_new": bool(b),
        "four_cells_comparable": False,
        "four_cells_equal": "not stored for nonestimable Stage6 row",
        "available_counts_equal": False,
        "estimability_equal": False,
        "BY_q_equal_within_rounding": False,
    }
    if a and b:
        columns = [
            ("a_target_case_reports", "a"),
            ("b_target_noncase_reports", "b"),
            ("c_background_case_reports", "c"),
            ("d_background_noncase_reports", "d"),
        ]
        r["four_cells_comparable"] = all(a.get(x) != "" for x, y in columns)
        if r["four_cells_comparable"]:
            r["four_cells_equal"] = all(int(a[x]) == int(b[y]) for x, y in columns)
        r["available_counts_equal"] = int(a["cohort_reports_n"]) == int(b["cohort_n"]) and int(
            a["target_reports_n"]
        ) == int(b["target_n"])
        r["estimability_equal"] = a["estimable"] == b["estimable"]
        if a["estimable"] == "True" and b["estimable"] == "True":
            r["old_BY_q"] = a["BY_adjusted_q_within_family"]
            r["new_BY_q"] = b["BY_q"]
            r["BY_q_equal_within_rounding"] = math.isclose(
                float(r["old_BY_q"]), float(r["new_BY_q"]), rel_tol=1e-6, abs_tol=1e-250
            )
        else:
            r["BY_q_equal_within_rounding"] = True
    r["pass"] = (
        (not r["four_cells_comparable"] or r["four_cells_equal"])
        and r["available_counts_equal"]
        and r["estimability_equal"]
        and r["BY_q_equal_within_rounding"]
    )
    rows.append(r)
write("annual_stage6_stage8_reconciliation.csv", rows)

oldcomp = {r["drug"]: r for r in read(STAGE6 / "阶段六_复合结局分解.csv") if r["database"] == "FAERS"}
newcomp = read(OUT / "broad_domain_contribution.csv")
ex = read(OUT / "broad_15_exclusive_combinations.csv")
inter = read(OUT / "broad_all_domain_intersections.csv")
comprows = []
for d, oldrow in oldcomp.items():
    nc = {r["domain"]: int(r["domain_n"]) for r in newcomp if r["drug"] == d}
    overlap = int(
        next(
            r["n"]
            for r in inter
            if r["drug"] == d and r["intersection"] == "dehydration & electrolyte"
        )
    )
    calc = nc["dehydration"] + nc["electrolyte"] - overlap
    checks = {
        "primary": (nc["primary"], int(oldrow["primary_union_reports_n"])),
        "renal": (nc["renal"], int(oldrow["renal_injury_domain_reports_n"])),
        "dehydration_electrolyte_union": (
            calc,
            int(oldrow["dehydration_electrolyte_domain_reports_n"]),
        ),
        "broad_union": (
            sum(int(r["n"]) for r in ex if r["drug"] == d),
            int(oldrow["broad_composite_reports_n"]),
        ),
    }
    for k, (obs, expected) in checks.items():
        comprows.append(
            {
                "drug": d,
                "check": k,
                "observed": obs,
                "original_stage6": expected,
                "pass": obs == expected,
            }
        )
write("composite_stage6_stage8_reconciliation.csv", comprows)

# Independently reconstruct BY with sorted ranks and reverse cumulative minima.
sens = read(OUT / "sensitivity_exact_ROR_BY.csv")
allrows = sens + newtime
qchecks = []
for family in sorted({r["family"] for r in allrows}):
    group = [r for r in allrows if r["family"] == family and r["estimable"] == "True"]
    m = len(group)
    p = np.array([float(r["Fisher_P"]) for r in group])
    order = np.argsort(p, kind="stable")
    h = sum(1 / k for k in range(1, m + 1))
    qs = np.minimum.accumulate((p[order] * m * h / np.arange(1, m + 1))[::-1])[::-1].clip(0, 1)
    actual = np.array([float(group[i]["BY_q"]) for i in order])
    err = float(np.max(np.abs(qs - actual))) if m else 0
    qchecks.append(
        {
            "family": family,
            "estimable_tests": m,
            "maximum_absolute_BY_difference": err,
            "pass": bool(np.allclose(qs, actual, rtol=1e-12, atol=1e-250)),
        }
    )
write("independent_BY_verification.csv", qchecks)

yearlist = list(newyears.values())
estim = [r for r in yearlist if r["estimable"] == "True"]
sig = [r for r in estim if float(r["BY_q"]) < 0.05]
summary = {
    "annual_possible_tests": len(yearlist),
    "annual_estimable_tests": len(estim),
    "annual_estimable_2019_onwards": sum(int(r["time_value"]) >= 2019 for r in estim),
    "annual_BY_below_0_05": len(sig),
    "annual_significant_higher_ROR": sum(float(r["ROR"]) > 1 for r in sig),
    "annual_significant_lower_ROR": sum(float(r["ROR"]) < 1 for r in sig),
    "annual_rows_with_original_four_cells": sum(r["four_cells_comparable"] for r in rows),
    "annual_all_available_counts_match_stage6": all(r["available_counts_equal"] for r in rows),
    "annual_all_four_cells_match_stage6": all(
        r["four_cells_equal"] for r in rows if r["four_cells_comparable"]
    ),
    "annual_all_BY_match_stage6": all(r["BY_q_equal_within_rounding"] for r in rows),
    "composite_all_match_stage6": all(r["pass"] for r in comprows),
    "independent_BY_all_pass": all(r["pass"] for r in qchecks),
}
(OUT / "independent_validation_summary.json").write_text(
    json.dumps(summary, indent=2), encoding="utf-8"
)
print(json.dumps(summary, indent=2))
if not all(
    summary[k]
    for k in [
        "annual_all_four_cells_match_stage6",
        "annual_all_BY_match_stage6",
        "composite_all_match_stage6",
        "independent_BY_all_pass",
    ]
):
    raise SystemExit("Validation difference found: inspect before integration")
