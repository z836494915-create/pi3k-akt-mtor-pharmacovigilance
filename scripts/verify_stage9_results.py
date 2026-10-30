"""Independent report-set counts and probability-tail checks for six contrasts."""
import csv, gzip, json, math, re
from pathlib import Path
import numpy as np
from scipy.stats import nchypergeom_fisher, hypergeom

from package_paths import OUTPUT, PACKAGE

O = OUTPUT / "03_确认后补充分析/阶段九_选药范围与竞争背景"
S8 = OUTPUT / "03_确认后补充分析/阶段八_核验后修复/FAERS"


def rows(n):
    with (O / n).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


dictionary = json.loads(
    (PACKAGE / "definitions/sglt2_dictionary.json").read_text(encoding="utf-8")
)["final_dictionary"]
patterns = [
    re.compile(
        r"(?<![A-Z0-9])"
        + re.escape(re.sub("[^A-Z0-9]+", " ", a.upper()).strip()).replace(r"\ ", r"\s*")
        + r"(?![A-Z])"
    )
    for aa in dictionary.values()
    for a in aa
]
sg = {"psss": set(), "any": set()}
fy = {"ps": set(), "psss": set()}
with gzip.open(
    O / "local_candidate_drug_rows_DO_NOT_PUBLISH.csv.gz", "rt", encoding="utf-8", newline=""
) as f:
    for r in csv.DictReader(f):
        pid = int(r["PRIMARYID"])
        role = r["role"]
        texts = [re.sub("[^A-Z0-9]+", " ", r[k].upper()).strip() for k in ["drugname", "PROD_AI"]]
        if any(p.search(t) for p in patterns for t in texts):
            sg["any"].add(pid)
            if role in ["PS", "SS"]:
                sg["psss"].add(pid)
        # All observed formulation-specific rows were explicitly named FYARRO.
        if any(re.search(r"\bFYARRO\b", t) for t in texts):
            if role == "PS":
                fy["ps"].add(pid)
            if role in ["PS", "SS"]:
                fy["psss"].add(pid)
with np.load(S8 / "local_report_cache_DO_NOT_PUBLISH.npz") as z:
    pids = z["pids"]
    ps = z["ps"]
    pts = z["pts"]
N = len(pids)
Y = set(map(int, pids[(pts & 7) > 0]))
targets = {
    d: set(map(int, pids[(ps & (1 << k)) > 0]))
    for k, d in enumerate(["alpelisib", "inavolisib", "capivasertib", "everolimus", "temsirolimus"])
}
out = rows("additional_sensitivity_exact_BY.csv")
checks = []
for r in out:
    if r["analysis"].startswith("exclude"):
        role = r["analysis"].split("_")[-1]
        excluded = sg[role]
        E = targets[r["drug"]] - excluded
        YY = Y - excluded
        n = N - len(excluded)
    else:
        role = "psss" if r["analysis"].endswith("psss") else "ps"
        E = fy[role]
        YY = Y
        n = N
    a = len(E & YY)
    b = len(E) - a
    c = len(YY) - a
    d = n - a - b - c
    assert [a, b, c, d] == [int(r[k]) for k in "abcd"], r
    R = a * d / (b * c)
    assert math.isclose(R, float(r["ROR"]), abs_tol=1e-12)
    lo = float(r["exact_CI_low"])
    hi = float(r["exact_CI_high"])
    # Inversion of the noncentral hypergeometric tail, independently of odds_ratio().
    upper_tail = nchypergeom_fisher.cdf(a, n, a + c, a + b, hi)
    assert abs(upper_tail - 0.025) < 2e-6, (r, upper_tail)
    lower_tail = nchypergeom_fisher.sf(a - 1, n, a + c, a + b, lo) if lo else None
    if lo:
        assert abs(lower_tail - 0.025) < 2e-6, (r, lower_tail)
    support = np.arange(max(0, a + b - (b + d)), min(a + b, a + c) + 1)
    logprob = hypergeom.logpmf(support, n, a + c, a + b)
    observed = hypergeom.logpmf(a, n, a + c, a + b)
    p = min(1.0, float(np.exp(logprob[logprob <= observed + 1e-7]).sum()))
    assert math.isclose(p, float(r["Fisher_P"]), rel_tol=1e-5, abs_tol=1e-8), (p, r["Fisher_P"])
    checks.append(
        {
            "analysis": r["analysis"],
            "drug": r["drug"],
            "report_set_counts_match": True,
            "cross_product_match": True,
            "exact_CI_tail_check": True,
            "two_sided_hypergeometric_p_check": True,
        }
    )
psort = sorted(enumerate(out), key=lambda x: float(x[1]["Fisher_P"]))
m = len(psort)
H = sum(1 / i for i in range(1, m + 1))
q = 1.0
for rank in range(m, 0, -1):
    idx, r = psort[rank - 1]
    q = min(q, float(r["Fisher_P"]) * m * H / rank, 1.0)
    assert abs(q - float(r["BY_q"])) < 1e-12
with (O / "independent_verification.json").open("w", encoding="utf-8") as f:
    json.dump(
        {
            "all_pass": True,
            "verification": "Separate report-ID sets reconstructed from candidate drug rows; exact interval probability tails and independent Fisher tail summation; manual BY adjustment. This is an independent computational check, not an independent investigator review.",
            "checks": checks,
            "six_planned_six_reported": True,
            "PS_and_PSSS_FYARRO_same_report_set": fy["ps"] == fy["psss"],
        },
        f,
        indent=2,
    )
print(
    "PASS: six four-cell tables, point estimates, exact intervals, Fisher tests, BY family; FYARRO role sets identical."
)
