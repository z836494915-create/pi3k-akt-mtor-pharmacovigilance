"""Numerical density integration independently checks Gamma quantile outputs."""
import argparse, csv, json, math, hashlib
from pathlib import Path
from collections import defaultdict
from scipy.integrate import quad
from scipy.special import gammaln


def rows(p):
    with p.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def posterior_cdf(x, h, r):
    # Log-coordinate integral, centered/scaled to avoid high-shape underflow.
    sq = math.sqrt(h)
    upper = sq * (math.log(x) - math.log(h / r))

    def density(z):
        u = z / sq
        if u > 700:
            return 0.0
        logp = h * (math.log(h) + u) - h * math.exp(u) - gammaln(h) - math.log(sq)
        return math.exp(logp) if logp > -745 else 0.0

    return quad(density, -math.inf, upper, epsabs=2e-10, epsrel=2e-10, limit=250)[0]


def verify(folder):
    data = rows(folder / "gamma_poisson_all_priors.csv")
    base = rows(folder / "gamma_poisson_default_prior.csv")
    assert len(data) == 513 and len(base) == 171
    errors = []
    groups = defaultdict(list)
    tested = 0
    na = 0
    closed = 0
    for x in data:
        key = tuple(x[k] for k in ["family", "database", "drug", "outcome", "analysis"])
        groups[key].append(x)
        a, b, c, d = [int(x[k]) for k in "abcd"]
        N = a + b + c + d
        E = (a + b) * (a + c) / N
        k = float(x["prior_k"])
        assert N == int(x["total_eligible_reports_n"]) and c + d == int(x["comparator_reports_n"])
        assert abs(E - float(x["expected_count"])) < 1e-10
        if x["estimable"] == "False":
            na += 1
            assert all(
                x[v] == ""
                for v in ["shrunken_OE", "CrI95_low", "CrI95_high", "nominal_lower_bound_gt_1"]
            )
            assert (a + b) == 0
            continue
        h = a + k
        r = E + k
        point = h / r
        assert math.isclose(point, float(x["shrunken_OE"]), rel_tol=1e-13)
        assert math.isclose(
            math.log2(point), float(x["log2_shrunken_OE"]), rel_tol=1e-13, abs_tol=1e-13
        )
        oe = a / E
        assert min(1, oe) - 1e-12 <= point <= max(1, oe) + 1e-12
        for col, p in [("CrI95_low", 0.025), ("CrI95_high", 0.975)]:
            bound = float(x[col])
            error = abs(posterior_cdf(bound, h, r) - p)
            errors.append(error)
            assert error < 1e-7, (key, k, col, error)
            if a == 0 and k == 1:
                assert math.isclose(bound, -math.log1p(-p) / r, rel_tol=1e-12)
                closed += 1
        assert (float(x["CrI95_low"]) > 1) == (x["nominal_lower_bound_gt_1"] == "True")
        tested += 1
    assert len(groups) == 171 and tested == 480 and na == 33
    for key, x in groups.items():
        assert sorted(float(r["prior_k"]) for r in x) == [0.25, 0.5, 1.0]
    summaries = rows(folder / "posterior_frequentist_descriptive_agreement.csv")
    for s in summaries:
        z = [
            x
            for x in data
            if x["family"] == s["family"]
            and float(x["prior_k"]) == float(s["prior_k"])
            and x["estimable"] == "True"
        ]
        p = [x["nominal_lower_bound_gt_1"] == "True" for x in z]
        f = [x["original_frequentist_elevated_screen"] == "True" for x in z]
        wanted = {
            "estimable_n": len(z),
            "nominal_posterior_lower_gt_1_n": sum(p),
            "frequentist_elevated_screen_n": sum(f),
            "both_n": sum(a and b for a, b in zip(p, f)),
            "posterior_only_n": sum(a and not b for a, b in zip(p, f)),
            "frequentist_only_n": sum(not a and b for a, b in zip(p, f)),
        }
        assert all(int(s[k]) == v for k, v in wanted.items())
    changes = [
        key for key, g in groups.items() if len(set(x["nominal_lower_bound_gt_1"] for x in g)) > 1
    ]
    discord = [
        {
            k: x[k]
            for k in [
                "family",
                "database",
                "drug",
                "outcome",
                "a",
                "reference_adjusted_P_or_q",
                "CrI95_low",
                "CrI95_high",
            ]
        }
        for x in base
        if x["estimable"] == "True"
        and x["nominal_lower_bound_gt_1"] != x["original_frequentist_elevated_screen"]
    ]
    assert len(changes) == 4 and len(discord) == 3
    result = {
        "all_pass": True,
        "scope": "Independent numerical method verification, not an independent analyst review or fresh raw-data extraction",
        "comparison_prior_rows": 513,
        "estimable_rows": tested,
        "nonestimable_rows_retained": na,
        "gamma_interval_endpoints_checked_by_density_integration": len(errors),
        "maximum_CDF_absolute_error": max(errors),
        "zero_count_exponential_closed_form_endpoints_checked": closed,
        "all_nine_descriptive_summaries_recomputed": True,
        "all_prior_boundary_changes": changes,
        "all_default_prior_frequentist_disagreements": discord,
        "inputs_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in folder.glob("*.csv")
        },
    }
    (folder / "independent_numerical_verification.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--results", type=Path, required=True)
    verify(p.parse_args().results)
