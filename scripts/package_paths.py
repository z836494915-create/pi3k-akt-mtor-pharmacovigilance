"""Single configuration entry point; no imports from the author's working tree."""
from __future__ import annotations
import json
import os
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1]
CONFIG_PATH = Path(os.environ.get("REVISION_CONFIG", str(PACKAGE / "config.local.json"))).resolve()
if not CONFIG_PATH.exists():
    raise RuntimeError(
        "Create config.local.json from config.example.json, or set REVISION_CONFIG to its absolute path"
    )
CONFIG = json.loads(CONFIG_PATH.read_text(encoding="utf-8-sig"))


def configured_path(key: str) -> Path:
    value = Path(CONFIG[key]).expanduser()
    return (CONFIG_PATH.parent / value).resolve() if not value.is_absolute() else value.resolve()


OUTPUT = configured_path("output_dir")
FAERS_ZIPS = configured_path("faers_zips")
FAERS_DB = configured_path("faers_db")
CANADA = configured_path("canada_dir")
JADER = configured_path("jader_dir")
STAGE3 = configured_path("audit_reference_dir") if CONFIG.get("audit_reference_dir") else OUTPUT / "optional_prior_audit"
