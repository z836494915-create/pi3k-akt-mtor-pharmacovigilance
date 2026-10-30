"""Verify the code distribution or run analyses from local source snapshots."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def check_code():
    entries = (ROOT / "manifest.sha256").read_text(encoding="utf-8").splitlines()
    seen = set()
    compiled = 0
    for entry in entries:
        expected, name = entry.split("  ", 1)
        path = (ROOT / name).resolve()
        if name in seen or not path.is_relative_to(ROOT) or not path.is_file():
            raise ValueError(f"Invalid distribution entry: {name}")
        seen.add(name)
        if digest(path) != expected:
            raise ValueError(f"File checksum differs: {name}")
        if path.suffix == ".py":
            compile(path.read_text(encoding="utf-8-sig"), str(path), "exec")
            compiled += 1
    return {"files_verified": len(entries), "python_files_compiled": compiled}


def read_config(path):
    data = json.loads(path.read_text(encoding="utf-8-sig"))
    paths = {}
    for name in ["faers_zips", "faers_db", "canada_dir", "jader_dir", "output_dir"]:
        value = Path(data[name]).expanduser()
        paths[name] = (path.parent / value).resolve() if not value.is_absolute() else value.resolve()
    return paths


def check_sources(paths):
    bases = {"FAERS": paths["faers_zips"], "Canada Vigilance": paths["canada_dir"],
             "JADER": paths["jader_dir"]}
    with (ROOT / "definitions/source_files.csv").open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError("Empty source-file manifest")
    checked = []
    for row in rows:
        path = bases[row["database"]] / row["filename"]
        if (not path.is_file() or path.stat().st_size != int(row["bytes"])
                or digest(path) != row["sha256"]):
            raise ValueError(f"Missing or changed source snapshot: {row['database']}/{row['filename']}")
        checked.append({"database": row["database"], "filename": row["filename"],
                        "sha256": row["sha256"]})
    return checked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["verify", "full"])
    parser.add_argument("--config", type=Path)
    args = parser.parse_args()
    if sys.flags.optimize:
        parser.error("Run without -O; the analytical scripts include assertions")
    integrity = check_code()
    if args.mode == "verify":
        if args.config:
            parser.error("verify does not use a local source configuration")
        print(json.dumps({"status": "passed", **integrity}, indent=2))
        return
    if args.config is None:
        parser.error("full requires --config")
    config = args.config.resolve()
    paths = read_config(config)
    output = paths["output_dir"]
    if output == ROOT or (output.exists() and any(output.iterdir())):
        raise FileExistsError("Choose a new or empty output directory")
    if paths["faers_db"].exists():
        raise FileExistsError("Choose a new FAERS database path")
    for source_dir in [paths["faers_zips"], paths["canada_dir"], paths["jader_dir"]]:
        if output.is_relative_to(source_dir) or paths["faers_db"].is_relative_to(source_dir):
            raise ValueError("Outputs must be outside the raw source directories")
    sources = check_sources(paths)
    output.mkdir(parents=True, exist_ok=True)
    logs = output / "execution_logs"
    logs.mkdir()
    env = dict(os.environ, REVISION_CONFIG=str(config), PYTHONDONTWRITEBYTECODE="1", PYTHONUTF8="1")
    env["PYTHONPATH"] = str(ROOT / "scripts") + os.pathsep + env.get("PYTHONPATH", "")
    status = {"status": "running", "source_checks": sources, "integrity": integrity, "steps": []}
    status_path = output / "analysis_status.json"

    def record():
        status_path.write_text(json.dumps(status, ensure_ascii=False, indent=2), encoding="utf-8")

    def run(script, *arguments):
        print("Running " + script, flush=True)
        start = time.monotonic()
        log_path = logs / (str(len(status["steps"])) + "_" + script.replace("/", "_") + ".log")
        with log_path.open("w", encoding="utf-8") as log:
            result = subprocess.run([sys.executable, str(ROOT / script), *map(str, arguments)],
                                    cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        status["steps"].append({"script": script, "returncode": result.returncode,
                                "seconds": round(time.monotonic()-start, 2)})
        record()
        if result.returncode:
            raise RuntimeError(f"Analysis failed: {script}; inspect execution_logs")

    record()
    try:
        for script in ["build_faers_demo.py", "run_stage5_primary_analysis.py",
                       "run_stage6_sensitivity.py", "run_stage6_reviewer_addendum.py",
                       "verify_and_summarize_stage6.py", "run_faers_stage8.py",
                       "verify_stage8_results.py", "run_external_revision_descriptives.py",
                       "run_background_analysis.py"]:
            run("scripts/" + script)
        primary = output / "03_确认后补充分析/阶段五_主要分析"
        secondary = output / "03_确认后补充分析/阶段六_敏感性分析"
        background = output / "03_确认后补充分析/阶段九_选药范围与竞争背景"
        shrinkage = output / "shrinkage"
        shrinkage_args = ["--primary", primary / "阶段五_主要结局ROR_Fisher_Holm.csv",
                          "--secondary", primary / "阶段五_次要结局ROR_Fisher_BY.csv",
                          "--additional", background / "additional_sensitivity_exact_BY.csv",
                          "--renal-reference", secondary / "阶段六_审稿人补充分析_renal_PT_ROR.csv",
                          "--output", shrinkage]
        run("scripts/run_shrinkage.py", "audit", *shrinkage_args)
        run("scripts/run_shrinkage.py", "analyze", *shrinkage_args)
        run("scripts/verify_shrinkage.py", "--results", shrinkage)
        jader = output / "clinical/JADER"
        canada = output / "clinical/Canada"
        shared = ROOT / "scripts/run_stage5_primary_analysis.py"
        run("clinical/JADER/audit_jader_clinical_availability.py", "--source-dir", paths["jader_dir"],
            "--stage5-script", shared, "--cohort-reference",
            output / "03_确认后补充分析/阶段八_核验后修复/External/external_stage5_reconciliation.csv",
            "--out", jader)
        run("clinical/JADER/describe_jader_PT_onset.py", "--source-dir", paths["jader_dir"],
            "--stage5-script", shared, "--audit-dir", jader, "--out", jader)
        shutil.copyfile(shared, primary / shared.name)
        run("clinical/Canada/audit_canada_clinical_fields.py", "--canada", paths["canada_dir"],
            "--stage5", primary, "--out", canada)
        run("clinical/Canada/run_canada_duplicate_flag_sensitivity.py", "--canada", paths["canada_dir"],
            "--stage5", primary, "--out", canada)
        status["status"] = "passed"
    except Exception:
        status["status"] = "failed"
        record()
        raise
    record()


if __name__ == "__main__":
    main()
