# Pharmacovigilance data processing and statistical analysis

Python scripts for FAERS, Canada Vigilance and JADER report-level data processing and statistical analysis. The repository contains executable code, drug/outcome dictionaries, a source-file checksum manifest and configuration instructions. It does not distribute source reports, precomputed analysis tables, figures, manuscripts or narrative-result generation code.

## Run

Use Python 3.12.3 and install the pinned numerical dependencies:

```text
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements-lock.txt
.venv\Scripts\python run_analysis.py verify
```

On macOS/Linux use `.venv/bin/python`; execution has been checked on Windows only. Run without `-O` because numerical checks include assertions.

Obtain the source snapshots independently: FAERS 2019Q1–2026Q1, Canada Vigilance `cvponline_extract_20260228` and JADER `202605`. The clinical-field scripts also require Canadian `report_links.txt` and JADER `hist202605.csv`. Exact filenames, sizes and hashes are in `definitions/source_files.csv`. There is no automatic download.

Copy `config.example.json` to `config.local.json`, enter local input paths, a new FAERS database path and an empty output directory, then run:

```text
.venv\Scripts\python run_analysis.py full --config config.local.json
```

The runner checks all source fingerprints before processing. It builds the FAERS latest-case table; applies drug and outcome matching; calculates primary, secondary, sensitivity and background comparisons; computes fixed-prior Gamma–Poisson estimates; and runs the Canadian/JADER clinical-field analyses. Each downstream analysis reads counts rebuilt in the same run. Numerical outputs and report-level caches remain local under the configured output directory. Stop and inspect `analysis_status.json` and `execution_logs/` if a step fails.

## Definitions

- FAERS uses primary-suspect coding for the primary analysis and PS/SS coding for sensitivity analyses; the other databases use native suspect roles.
- Main estimates use report-level four-cell tables, cross-product ROR, conditional exact confidence intervals and two-sided Fisher tests. Primary comparisons use Holm adjustment; exploratory families use Benjamini–Yekutieli adjustment.
- Drug names, JADER preferred-term mappings and outcome unions are supplied in `definitions/`. Background exclusions use explicit ingredient/formulation dictionaries. Missing information remains unknown.
- Fixed-prior Gamma–Poisson calculations use shape/rate values of 0.25, 0.5 and 1; they are complementary to the frequentist calculations.
- Prior cohort-flow summaries are optional local audit references (`audit_reference_dir` in the configuration). When absent, the primary script records that the historical comparison was not performed. Report sets and statistical counts are rebuilt from the specified raw sources, and subsequent stages reconcile against that run.

The scripts are tied to the listed snapshots and include numerical consistency checks for those releases. Adapting them to new source versions requires a new audit. These are report-level analyses; they do not establish patient incidence, causal effects or verified clinical diagnoses. Computational checks do not constitute clinical adjudication.

## Availability and rights

Public release: `v1.0.0`. No software license has been granted. Public visibility does not itself grant a general license to reuse, modify or redistribute the code. Source-database terms and third-party dependencies retain their own terms.
