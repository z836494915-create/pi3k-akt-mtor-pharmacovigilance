from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path

from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
from statsmodels.stats.multitest import multipletests


from package_paths import OUTPUT

OUT = OUTPUT / "03_确认后补充分析" / "阶段六_敏感性分析"
STAGE5 = OUT.parent / "阶段五_主要分析"


def read_rows(name, root=OUT):
    with (root / name).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_rows(name, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def near(a, b, rel=2e-4, abs_=1e-10):
    return math.isclose(float(a), float(b), rel_tol=rel, abs_tol=abs_)


def sig(x):
    try:
        return f"{float(x):.3g}"
    except (TypeError, ValueError):
        return "—"


def main():
    sens = read_rows("阶段六_敏感性ROR结果.csv")
    checks = []
    estimable = [r for r in sens if r.get("estimable") == "True"]
    cell_errors, ror_errors, p_errors, ci_errors = [], [], [], []
    p_abs_max = 0.0
    for r in estimable:
        a, b, c, d = [
            int(float(r[k]))
            for k in (
                "a_target_case_reports",
                "b_target_noncase_reports",
                "c_background_case_reports",
                "d_background_noncase_reports",
            )
        ]
        n = int(float(r["cohort_reports_n"]))
        target_n = int(float(r["target_reports_n"]))
        if a + b + c + d != n or a + b != target_n:
            cell_errors.append((r["sensitivity"], r["drug"], r.get("report_year_by_FDA_DT")))
        cross = (a * d) / (b * c) if b * c else (math.inf if a * d else math.nan)
        saved_ror = r.get("ROR_cross_product", "")
        if not saved_ror or not math.isclose(cross, float(saved_ror), rel_tol=2e-4, abs_tol=1e-10):
            ror_errors.append((r["sensitivity"], r["drug"], r.get("report_year_by_FDA_DT")))
        p = float(fisher_exact([[a, b], [c, d]], alternative="two-sided").pvalue)
        saved_p = float(r["fisher_two_sided_p"])
        p_abs_max = max(p_abs_max, abs(p - saved_p))
        if not math.isclose(p, saved_p, rel_tol=2e-4, abs_tol=1e-5):
            p_errors.append(
                (r["sensitivity"], r["drug"], r.get("report_year_by_FDA_DT"), p, saved_p)
            )
        ci = odds_ratio([[a, b], [c, d]], kind="conditional").confidence_interval(
            confidence_level=0.95
        )
        lo, hi = float(ci.low), float(ci.high)
        if not (
            near(lo, r["exact_conditional_95CI_lower"])
            and near(hi, r["exact_conditional_95CI_upper"])
        ):
            ci_errors.append((r["sensitivity"], r["drug"], r.get("report_year_by_FDA_DT")))
    checks.extend(
        [
            {
                "check": "All estimable 2x2 tables sum to cohort and target totals",
                "n_checked": len(estimable),
                "n_failed": len(cell_errors),
                "passed": not cell_errors,
                "details": str(cell_errors[:5]),
            },
            {
                "check": "Cross-product ROR independently recalculated",
                "n_checked": len(estimable),
                "n_failed": len(ror_errors),
                "passed": not ror_errors,
                "details": str(ror_errors[:5]),
            },
            {
                "check": "Two-sided Fisher P independently recalculated (allowing saved 5-significant-digit rounding)",
                "n_checked": len(estimable),
                "n_failed": len(p_errors),
                "passed": not p_errors,
                "details": f"max absolute saved-rounding difference={p_abs_max:.3g}; failures={p_errors[:5]}",
            },
            {
                "check": "Conditional exact 95% CI independently recalculated",
                "n_checked": len(estimable),
                "n_failed": len(ci_errors),
                "passed": not ci_errors,
                "details": str(ci_errors[:5]),
            },
        ]
    )

    by_family = defaultdict(list)
    for r in estimable:
        by_family[r["sensitivity_family"]].append(r)
    q_errors = []
    for family, group in by_family.items():
        qs = multipletests([float(r["fisher_two_sided_p"]) for r in group], method="fdr_by")[1]
        for r, q in zip(group, qs):
            if not near(q, r["BY_adjusted_q_within_family"], rel=8e-4, abs_=1e-9):
                q_errors.append(
                    (
                        family,
                        r["drug"],
                        r.get("report_year_by_FDA_DT"),
                        q,
                        r["BY_adjusted_q_within_family"],
                    )
                )
    checks.append(
        {
            "check": "BY q-values independently recalculated within sensitivity family",
            "n_checked": len(estimable),
            "n_failed": len(q_errors),
            "passed": not q_errors,
            "details": str(q_errors[:5]),
        }
    )

    # Confirm annual report universes aggregate exactly from mutually exclusive FDA_DT quarters.
    quarter_rows = read_rows("阶段六_FAERS季度报告计数.csv")
    quarter_n = {}
    for r in quarter_rows:
        quarter_n[r["report_quarter_by_FDA_DT"]] = int(r["quarter_reports_n"])
    annual = [r for r in sens if r["sensitivity"] == "FAERS_annual_by_FDA_DT"]
    annual_checks = []
    for r in annual:
        year = int(r["report_year_by_FDA_DT"])
        qsum = sum(n for q, n in quarter_n.items() if int(q[:4]) == year)
        if qsum != int(r["cohort_reports_n"]):
            annual_checks.append((year, r["drug"], qsum, r["cohort_reports_n"]))
    checks.append(
        {
            "check": "Annual report totals equal sum of FDA_DT quarter totals",
            "n_checked": len(annual),
            "n_failed": len(annual_checks),
            "passed": not annual_checks,
            "details": str(annual_checks[:5]),
        }
    )

    # Confirm primary report counts reproduce the locked Stage 5 analysis.
    stage5 = read_rows("阶段五_主要结局ROR_Fisher_Holm.csv", STAGE5)
    faers_sens_rows = {
        r["drug"]: r for r in sens if r["sensitivity"] == "FAERS_PS_vs_no_target_SS_only"
    }
    reconciliation = []
    for r in stage5:
        if r["database"] != "FAERS":
            continue
        s = faers_sens_rows[r["drug"]]
        matched = int(r["target_reports_n"]) == int(s["target_reports_n"]) and int(
            r["a_target_case_reports"]
        ) == int(s["a_target_case_reports"])
        reconciliation.append(matched)
    checks.append(
        {
            "check": "FAERS PS target/case counts match Stage 5",
            "n_checked": len(reconciliation),
            "n_failed": sum(not x for x in reconciliation),
            "passed": all(reconciliation),
            "details": "",
        }
    )

    spontaneous = read_rows("阶段六_Canada报告类型核验.csv")
    spont_n = sum(
        int(r["reports_n"]) for r in spontaneous if r["report_type"].casefold() == "spontaneous"
    )
    checks.append(
        {
            "check": "Canada spontaneous count matches Stage 3",
            "n_checked": spont_n,
            "n_failed": int(spont_n != 681308),
            "passed": spont_n == 681308,
            "details": "Stage 3 locked count=681308",
        }
    )
    write_rows("阶段六_独立复核.csv", checks)

    if not all(row["passed"] for row in checks):
        raise ValueError("Sensitivity verification failed")
    print(f"Checks passed: {len(checks)}")


if __name__ == "__main__":
    main()
