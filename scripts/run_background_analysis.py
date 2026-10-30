"""Audit drug-name candidates, then run the specified background comparisons."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import subprocess
import sys
from pathlib import Path

from package_paths import OUTPUT, PACKAGE

HERE = Path(__file__).resolve().parent
OUT = OUTPUT / "03_确认后补充分析/阶段九_选药范围与竞争背景"


def main():
    OUT.mkdir(parents=True, exist_ok=False)
    subprocess.run([sys.executable, str(HERE / "run_stage9_background.py"), "audit"], check=True)
    import run_stage9_background as background

    definition_path = PACKAGE / "definitions/sglt2_dictionary.json"
    definition = json.loads(definition_path.read_text(encoding="utf-8"))
    ambiguous = set(definition["ambiguous_names_excluded"])
    ambiguous_rows = 0
    conflicts = 0
    with gzip.open(OUT / "local_candidate_drug_rows_DO_NOT_PUBLISH.csv.gz", "rt",
                   encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            fields = [background.norm(row["drugname"]), background.norm(row["PROD_AI"])]
            hits = [{key for key, pattern in background.SG_PAT.items() if pattern.search(value)}
                    for value in fields]
            if hits[0] and hits[1] and not hits[0] & hits[1]:
                conflicts += 1
            if not (hits[0] | hits[1]) and "GLIFLOZIN" in " ".join(fields):
                if fields[0] not in ambiguous:
                    raise ValueError("An unreviewed ingredient name requires a new dictionary audit")
                ambiguous_rows += 1
    if conflicts:
        raise ValueError("Drug-name and active-ingredient classifications conflict")
    audit = json.loads((OUT / "audit_extraction.json").read_text(encoding="utf-8"))
    if not audit["unique_sorted_PRIMARYID"]:
        raise ValueError("Report index failed its uniqueness/order check")
    gate = {
        "analysis_allowed": True,
        "unresolved_candidate_decision_n": 0,
        "ambiguous_identity_excluded_rows_n": ambiguous_rows,
        "dictionary_amendment_sha256": hashlib.sha256(definition_path.read_bytes()).hexdigest(),
        "execution_mode": "new_DRUG_source_extraction",
    }
    (OUT / "audit_acceptance.json").write_text(json.dumps(gate, indent=2), encoding="utf-8")
    for name, arguments in [("run_stage9_background.py", ["analyze"]),
                            ("verify_stage9_results.py", [])]:
        subprocess.run([sys.executable, str(HERE / name), *arguments], check=True)


if __name__ == "__main__":
    main()
