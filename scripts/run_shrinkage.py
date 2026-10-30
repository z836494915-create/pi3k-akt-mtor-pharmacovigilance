"""Frozen-count Gamma-Poisson complement; no patient-level causal interpretation."""
from __future__ import annotations
import argparse, csv, hashlib, json, math
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
import scipy
from scipy.stats import gamma


def read(p):
    with Path(p).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write(p, rows):
    with Path(p).open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        w.writeheader()
        w.writerows(rows)


def save(p, obj):
    p.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding="utf-8")


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def source_paths(args):
    return [args.primary, args.secondary, args.additional, args.renal_reference]


def load_counts(args):
    allrows = []
    for family, path in [
        ("primary", args.primary),
        ("secondary", args.secondary),
        ("stage9_extra", args.additional),
    ]:
        for r in read(path):
            keys = (
                [
                    "a_target_case_reports",
                    "b_target_noncase_reports",
                    "c_background_case_reports",
                    "d_background_noncase_reports",
                ]
                if family != "stage9_extra"
                else list("abcd")
            )
            vals = [int(r[k]) for k in keys]
            assert all(str(v) == r[k].strip() and v >= 0 for v, k in zip(vals, keys))
            a, b, c, d = vals
            N = sum(vals)
            T = a + b
            Y = a + c
            key = (
                r.get("database", "FAERS"),
                r["drug"],
                r.get("outcome", "primary_union"),
                r.get("analysis", "original"),
            )
            qkey = (
                "holm_adjusted_p"
                if family == "primary"
                else "BY_adjusted_q"
                if family == "secondary"
                else "BY_q"
            )
            rorfield = "ROR_cross_product" if family != "stage9_extra" else "ROR"
            estimable = T > 0 and c + d > 0 and Y > 0 and b + d > 0
            assert estimable == (r["estimable"].lower() == "true")
            allrows.append(
                dict(
                    family=family,
                    database=key[0],
                    drug=key[1],
                    outcome=key[2],
                    analysis=key[3],
                    a=a,
                    b=b,
                    c=c,
                    d=d,
                    total_eligible_reports_n=N,
                    target_reports_n=T,
                    comparator_reports_n=c + d,
                    total_outcome_reports_n=Y,
                    estimable=estimable,
                    reference_ROR=r[rorfield],
                    reference_adjusted_P_or_q=r.get(qkey, ""),
                    reference_CI_low=r.get(
                        "exact_conditional_95CI_lower", r.get("exact_CI_low", "")
                    ),
                )
            )
    assert (
        len(allrows) == 171
        and len(
            {(r["family"], r["database"], r["drug"], r["outcome"], r["analysis"]) for r in allrows}
        )
        == 171
    )
    assert Counter(r["family"] for r in allrows) == {
        "primary": 15,
        "secondary": 150,
        "stage9_extra": 6,
    }
    return allrows


def audit(args):
    o = args.output
    o.mkdir(parents=True, exist_ok=True)
    assert not (o / "analysis_completion.json").exists(), "Use a fresh output directory."
    rr = load_counts(args)
    # The 45 renal rows are already present in the 150-row secondary table. Check, never append them.
    renal = read(args.renal_reference)
    idx = {(r["database"], r["drug"], r["outcome"]): r for r in rr if r["family"] == "secondary"}
    assert len(renal) == 45
    for r in renal:
        x = idx[r["database"], r["drug"], r["outcome"]]
        assert [x[k] for k in "abcd"] == [
            int(r[k])
            for k in [
                "a_target_case_reports",
                "b_target_noncase_reports",
                "c_background_case_reports",
                "d_background_noncase_reports",
            ]
        ]
    margins = defaultdict(set)
    pop = defaultdict(set)
    target = defaultdict(set)
    for r in rr:
        if r["family"] == "stage9_extra":
            continue
        margins[r["database"], r["outcome"]].add(r["total_outcome_reports_n"])
        pop[r["database"]].add(r["total_eligible_reports_n"])
        target[r["database"], r["drug"]].add(r["target_reports_n"])
    assert all(len(v) == 1 for d in [margins, pop, target] for v in d.values())
    # Legacy label is N, not the non-target count c+d. Preserve frozen source files and document the field.
    legacy = []
    for p in [args.primary, args.secondary, args.renal_reference]:
        for r in read(p):
            N = sum(
                int(r[k])
                for k in [
                    "a_target_case_reports",
                    "b_target_noncase_reports",
                    "c_background_case_reports",
                    "d_background_noncase_reports",
                ]
            )
            assert int(r["background_reports_n"]) == N
        legacy.append({"file": p.name, "legacy_background_reports_n_equals_total_N": True})
    sources = [
        {"file": p.name, "bytes": p.stat().st_size, "sha256": sha(p)} for p in source_paths(args)
    ]
    write(o / "audited_four_cell_inputs.csv", rr)
    plan = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "phase": "After four-cell audit and before new posterior estimates",
        "design": "Report-level case-noncase complement using the same frozen data; no new cohort, exposure or outcome definition.",
        "universe": "15 original primary + 150 original secondary + 6 stage-nine contrasts; all planned rows retained. The 45 renal rows already in secondary are checked but not appended.",
        "O": "a",
        "N": "a+b+c+d",
        "E": "(a+b)*(a+c)/N; entire database independence baseline including target reports",
        "likelihood": "O | lambda,E ~ Poisson(E*lambda), conditional on plug-in E",
        "prior": "lambda ~ Gamma(k,k), shape-rate parametrization, mean 1; not fitted to current results",
        "posterior": "Gamma(a+k,E+k)",
        "priors": [0.25, 0.5, 1.0],
        "default_prior": 0.5,
        "point": "posterior mean lambda=(a+k)/(E+k); log2 of this mean is the shrinkage log2 OE score, not the posterior mean of log2(lambda)",
        "interval": "Equal-tailed 95% posterior credible interval from numerical Gamma quantiles; log2 endpoints also supplied; no approximate normal/shrinkage limit formula.",
        "comparisons": "Show unshrunk O/E separately from ROR. Never label the difference between ROR and shrunken O/E as shrinkage alone. Do not label this BCPNN or fitted empirical-Bayes MGPS.",
        "nominal_flag": "Lower credible bound >1, described only as nominal posterior compatibility; no replacement of the original count/exact-CI/Holm or BY decisions.",
        "multiplicity": "Original families remain unchanged; these fixed-prior intervals are not multiplicity-adjusted and are not a new multiplicity solution.",
        "all_results": "All 171 rows x 3 fixed priors; nonestimable rows stay NA. No selective reporting of supportive combinations.",
        "scope_limits": [
            "Expected counts are treated as fixed; dependence, reporting bias and incomplete clinical data remain.",
            "OE and ROR differ as estimands, especially when marginal event reporting is common.",
            "Fixed-prior shrinkage can materially influence rare reports; test half/double prior strength.",
            "Neither a low value nor a zero event count establishes clinical protection or safety.",
        ],
        "sources": [
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC6331976/",
            "https://doi.org/10.1177/0962280211403604",
            "https://search.r-project.org/CRAN/refmans/pvda/html/ic.html",
        ],
        "input_manifest": sources,
    }
    p = o / "analysis_plan_before_posteriors.json"
    assert not p.exists(), "Do not overwrite a frozen plan."
    save(p, plan)
    gate = {
        "completed_at": datetime.now().isoformat(timespec="seconds"),
        "audit_pass": True,
        "rows": 171,
        "estimable": sum(r["estimable"] for r in rr),
        "unique_keys": True,
        "missing_count_cells": 0,
        "negative_count_cells": 0,
        "all_margins_consistent": True,
        "renal_duplicate_rows_checked_not_reappended": 45,
        "database_report_totals": {k: next(iter(v)) for k, v in pop.items()},
        "legacy_field_check": legacy,
        "sources": sources,
        "input_csv_sha256": sha(o / "audited_four_cell_inputs.csv"),
        "plan_sha256": sha(p),
        "analysis_allowed": True,
        "stage_transition": "Continuation and feasible supplementary analyses were authorized by the user; no new clinical model or public publication is undertaken.",
    }
    save(o / "data_audit_and_gate.json", gate)
    print(json.dumps(gate, ensure_ascii=False, indent=2))


def analyze(args):
    o = args.output
    gate = json.loads((o / "data_audit_and_gate.json").read_text(encoding="utf-8"))
    plan = json.loads((o / "analysis_plan_before_posteriors.json").read_text(encoding="utf-8"))
    assert (
        gate["analysis_allowed"]
        and gate["audit_pass"]
        and sha(o / "audited_four_cell_inputs.csv") == gate["input_csv_sha256"]
    )
    assert (
        sha(o / "analysis_plan_before_posteriors.json") == gate["plan_sha256"]
    ), "The frozen model/priors must not change after audit."
    for p, r in zip(source_paths(args), gate["sources"]):
        assert sha(p) == r["sha256"]
    rr = load_counts(args)
    output = []
    for r in rr:
        E = (
            r["target_reports_n"] * r["total_outcome_reports_n"] / r["total_eligible_reports_n"]
            if r["total_eligible_reports_n"]
            else 0.0
        )
        for k in plan["priors"]:
            row = {**r, "prior_k": k, "expected_count": E}
            if r["estimable"]:
                shape = r["a"] + k
                rate = E + k
                lo, hi = gamma.ppf([0.025, 0.975], a=shape, scale=1 / rate)
                oe = r["a"] / E
                mean = shape / rate
                logscore = math.log2(mean)
                row.update(
                    unshrunk_OE=oe,
                    posterior_shape=shape,
                    posterior_rate=rate,
                    shrunken_OE=mean,
                    CrI95_low=lo,
                    CrI95_high=hi,
                    log2_shrunken_OE=logscore,
                    log2_CrI95_low=math.log2(lo),
                    log2_CrI95_high=math.log2(hi),
                    posterior_probability_lambda_gt_1=float(gamma.sf(1, a=shape, scale=1 / rate)),
                    nominal_lower_bound_gt_1=bool(lo > 1),
                    shrunken_over_unshrunk=mean / oe if oe else "",
                    original_frequentist_elevated_screen=(
                        (r["family"] != "primary" or r["a"] >= 3)
                        and float(r["reference_CI_low"]) > 1
                        and float(r["reference_ROR"]) > 1
                        and float(r["reference_adjusted_P_or_q"]) < 0.05
                    ),
                )
                assert min(oe, 1) - 1e-12 <= mean <= max(oe, 1) + 1e-12
            else:
                row.update(
                    {
                        k: ""
                        for k in [
                            "unshrunk_OE",
                            "posterior_shape",
                            "posterior_rate",
                            "shrunken_OE",
                            "CrI95_low",
                            "CrI95_high",
                            "log2_shrunken_OE",
                            "log2_CrI95_low",
                            "log2_CrI95_high",
                            "posterior_probability_lambda_gt_1",
                            "nominal_lower_bound_gt_1",
                            "shrunken_over_unshrunk",
                            "original_frequentist_elevated_screen",
                        ]
                    }
                )
            output.append(row)
    write(o / "gamma_poisson_all_priors.csv", output)
    standard = [r for r in output if r["prior_k"] == 0.5]
    write(o / "gamma_poisson_default_prior.csv", standard)
    summaries = []
    for family in ["primary", "secondary", "stage9_extra"]:
        for k in plan["priors"]:
            s = [
                r for r in output if r["family"] == family and r["prior_k"] == k and r["estimable"]
            ]
            summaries.append(
                dict(
                    family=family,
                    prior_k=k,
                    estimable_n=len(s),
                    nominal_posterior_lower_gt_1_n=sum(r["nominal_lower_bound_gt_1"] for r in s),
                    frequentist_elevated_screen_n=sum(
                        r["original_frequentist_elevated_screen"] for r in s
                    ),
                    both_n=sum(
                        r["nominal_lower_bound_gt_1"] and r["original_frequentist_elevated_screen"]
                        for r in s
                    ),
                    posterior_only_n=sum(
                        r["nominal_lower_bound_gt_1"]
                        and not r["original_frequentist_elevated_screen"]
                        for r in s
                    ),
                    frequentist_only_n=sum(
                        not r["nominal_lower_bound_gt_1"]
                        and r["original_frequentist_elevated_screen"]
                        for r in s
                    ),
                )
            )
    write(o / "posterior_frequentist_descriptive_agreement.csv", summaries)
    sensitivity = []
    for r in standard:
        s = [
            x
            for x in output
            if (x["family"], x["database"], x["drug"], x["outcome"], x["analysis"])
            == (r["family"], r["database"], r["drug"], r["outcome"], r["analysis"])
        ]
        item = {
            k: r[k]
            for k in [
                "family",
                "database",
                "drug",
                "outcome",
                "analysis",
                "a",
                "expected_count",
                "estimable",
            ]
        }
        for x in s:
            tag = str(x["prior_k"])
            item.update(
                {
                    f"{tag}_" + key: x[key]
                    for key in [
                        "shrunken_OE",
                        "CrI95_low",
                        "CrI95_high",
                        "nominal_lower_bound_gt_1",
                    ]
                }
            )
        item["nominal_flag_changes_across_priors"] = (
            len({x["nominal_lower_bound_gt_1"] for x in s}) > 1
        )
        sensitivity.append(item)
    write(o / "prior_sensitivity_by_comparison.csv", sensitivity)
    save(
        o / "analysis_completion.json",
        {
            "completed_at": datetime.now().isoformat(timespec="seconds"),
            "rows": len(output),
            "estimable_comparison_prior_rows": sum(r["estimable"] for r in output),
            "scipy": scipy.__version__,
            "original_tests_unchanged": True,
            "data_audit_sha256": sha(o / "data_audit_and_gate.json"),
            "plan_sha256": sha(o / "analysis_plan_before_posteriors.json"),
        },
    )
    print(json.dumps(summaries, indent=2))
    print("Primary default-prior estimates:")
    for r in standard:
        if r["family"] == "primary":
            print(
                {
                    k: r.get(k)
                    for k in [
                        "database",
                        "drug",
                        "a",
                        "expected_count",
                        "unshrunk_OE",
                        "shrunken_OE",
                        "CrI95_low",
                        "CrI95_high",
                        "nominal_lower_bound_gt_1",
                    ]
                }
            )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("action", choices=["audit", "analyze"])
    for name in ["primary", "secondary", "additional", "renal-reference", "output"]:
        p.add_argument("--" + name, type=Path, required=True)
    a = p.parse_args()
    audit(a) if a.action == "audit" else analyze(a)
