"""FAERS DEMO construction: exact original deletion/latest-case rules.
Only IO safety/configuration and progress messages differ from the recorded source.
No source database is modified; an existing output is refused.
"""
from __future__ import annotations
import csv, io, json, re, sqlite3, zipfile
from pathlib import Path
from package_paths import FAERS_ZIPS, FAERS_DB, OUTPUT

RAW_DIR = FAERS_ZIPS
DB_PATH = FAERS_DB
START_Q = "2019q1"
END_Q = "2026q1"


def quarter_from_name(path: Path) -> str:
    m = re.search(r"faers_ascii_(\d{4}q[1-4])\.zip$", path.name, re.I)
    return m.group(1).lower() if m else path.stem.lower()


def quarter_sort_key(q: str) -> tuple[int, int]:
    m = re.match(r"(\d{4})q([1-4])", q.lower())
    if not m:
        return (9999, 9)
    return (int(m.group(1)), int(m.group(2)))


def list_zips() -> list[Path]:
    zips = sorted(
        RAW_DIR.glob("faers_ascii_*.zip"), key=lambda p: quarter_sort_key(quarter_from_name(p))
    )
    selected = []
    for z in zips:
        q = quarter_from_name(z)
        if quarter_sort_key(START_Q) <= quarter_sort_key(q) <= quarter_sort_key(END_Q):
            selected.append(z)
    return selected


def find_member(zf: zipfile.ZipFile, prefix: str) -> str | None:
    prefix = prefix.upper()
    for name in zf.namelist():
        base = Path(name).name.upper()
        if base.startswith(prefix) and base.endswith(".TXT"):
            return name
    return None


def find_deleted_members(zf: zipfile.ZipFile) -> list[str]:
    out = []
    for name in zf.namelist():
        lower = name.lower()
        base = Path(name).name.upper()
        if lower.endswith(".txt") and (
            "deleted/" in lower or base.startswith(("DELETE", "ALLDELETED", "ADR"))
        ):
            out.append(name)
    return out


def iter_table(zf: zipfile.ZipFile, member: str):
    with zf.open(member) as raw:
        text = io.TextIOWrapper(raw, encoding="latin-1", errors="replace", newline="")
        reader = csv.DictReader(text, delimiter="$")
        for row in reader:
            yield {k.lower(): (v or "").strip() for k, v in row.items() if k is not None}


def parse_int(value: str, default: int = 0) -> int:
    try:
        return int(value)
    except Exception:
        return default


def read_deleted_caseids(zips: list[Path]) -> set[str]:
    deleted = set()
    for zp in zips:
        with zipfile.ZipFile(zp) as zf:
            for member in find_deleted_members(zf):
                with zf.open(member) as raw:
                    text = io.TextIOWrapper(raw, encoding="latin-1", errors="replace")
                    for line in text:
                        caseid = line.strip().split("$")[0].strip()
                        if caseid and caseid.lower() != "caseid":
                            deleted.add(caseid)
    return deleted


def init_db(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        PRAGMA journal_mode = WAL;
        PRAGMA synchronous = OFF;
        PRAGMA temp_store = MEMORY;
        DROP TABLE IF EXISTS demo_latest;
        CREATE TABLE demo_latest (
            caseid TEXT PRIMARY KEY,
            primaryid TEXT NOT NULL,
            caseversion INTEGER NOT NULL,
            fda_dt TEXT,
            fda_dt_sort INTEGER,
            event_dt TEXT,
            init_fda_dt TEXT,
            age TEXT,
            age_cod TEXT,
            sex TEXT,
            reporter_country TEXT,
            occr_country TEXT,
            occp_cod TEXT,
            rept_cod TEXT,
            source_quarter TEXT
        );
        """
    )


def build_demo_latest(zips: list[Path], deleted: set[str]) -> dict[str, int]:
    stats = {"raw_demo_rows": 0, "deleted_demo_rows_skipped": 0}
    if DB_PATH.exists():
        raise FileExistsError(f"Refusing to overwrite existing database: {DB_PATH}")
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(DB_PATH)
    init_db(conn)
    sql = """
    INSERT INTO demo_latest (
        caseid, primaryid, caseversion, fda_dt, fda_dt_sort, event_dt, init_fda_dt,
        age, age_cod, sex, reporter_country, occr_country, occp_cod, rept_cod, source_quarter
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ON CONFLICT(caseid) DO UPDATE SET
        primaryid=excluded.primaryid,
        caseversion=excluded.caseversion,
        fda_dt=excluded.fda_dt,
        fda_dt_sort=excluded.fda_dt_sort,
        event_dt=excluded.event_dt,
        init_fda_dt=excluded.init_fda_dt,
        age=excluded.age,
        age_cod=excluded.age_cod,
        sex=excluded.sex,
        reporter_country=excluded.reporter_country,
        occr_country=excluded.occr_country,
        occp_cod=excluded.occp_cod,
        rept_cod=excluded.rept_cod,
        source_quarter=excluded.source_quarter
    WHERE
        excluded.caseversion > demo_latest.caseversion OR
        (excluded.caseversion = demo_latest.caseversion AND excluded.fda_dt_sort >= demo_latest.fda_dt_sort)
    """
    for zp in zips:
        q = quarter_from_name(zp)
        print(f"Building DEMO {q}", flush=True)
        batch = []
        with zipfile.ZipFile(zp) as zf:
            member = find_member(zf, "DEMO")
            if not member:
                continue
            for row in iter_table(zf, member):
                stats["raw_demo_rows"] += 1
                caseid = row.get("caseid", "")
                if not caseid or caseid in deleted:
                    stats["deleted_demo_rows_skipped"] += 1
                    continue
                fda_dt = row.get("fda_dt", "")
                batch.append(
                    (
                        caseid,
                        row.get("primaryid", ""),
                        parse_int(row.get("caseversion", "")),
                        fda_dt,
                        parse_int(fda_dt),
                        row.get("event_dt", ""),
                        row.get("init_fda_dt", ""),
                        row.get("age", ""),
                        row.get("age_cod", ""),
                        row.get("sex", ""),
                        row.get("reporter_country", ""),
                        row.get("occr_country", ""),
                        row.get("occp_cod", ""),
                        row.get("rept_cod", ""),
                        q,
                    )
                )
                if len(batch) >= 50000:
                    conn.executemany(sql, batch)
                    conn.commit()
                    batch.clear()
        if batch:
            conn.executemany(sql, batch)
            conn.commit()
    conn.execute("CREATE INDEX idx_demo_primaryid ON demo_latest(primaryid)")
    conn.commit()
    stats["latest_reports_after_deleted"] = conn.execute(
        "SELECT COUNT(*) FROM demo_latest"
    ).fetchone()[0]
    conn.close()
    return stats


def main():
    zips = list_zips()
    if len(zips) != 29:
        raise RuntimeError(f"Expected 29 quarterly ZIPs; found {len(zips)}")
    deleted = read_deleted_caseids(zips)
    stats = build_demo_latest(zips, deleted)
    stats["deleted_caseids_in_source_lists"] = len(deleted)
    stats["quarter_files"] = len(zips)
    stats[
        "rule"
    ] = "Exclude all listed CASEIDs across selected releases; choose highest CASEVERSION, then greatest numeric FDA_DT; later encountered row wins a complete tie. Original quarterly order retained."
    if stats["latest_reports_after_deleted"] != 10463182:
        raise RuntimeError(f"Reconstructed count differs from locked release: {stats}")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / "faers_demo_construction.json").write_text(
        json.dumps(stats, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(stats, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
