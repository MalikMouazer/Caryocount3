#!/usr/bin/env python3
import argparse
import io
import math
import os
import shutil
import subprocess
import types
from datetime import datetime
import sys
from typing import Callable, Optional, Tuple

import pandas as pd

from My_expert_karyo_functions import analyser_formule

DEFAULT_SHEET_URL = (
    "https://docs.google.com/spreadsheets/d/"
    "1_QPAzEr3PNHaNVu8Qxuv_hWCUoBQqqRcNSnX-Q0Egm8/edit?gid=98233212#gid=98233212"
)
DEFAULT_LOCAL_FILE = "comptage_local_MYC.xlsx"

# Baseline du fichier local:
# - changer ELISE_FILE_DATE quand Elise envoie un nouveau comptage manuel;
# - changer FROZEN_RULES_REF quand on décide de figer une nouvelle version
#   des règles. Utiliser un hash de commit Git, un tag, ou "HEAD".
ELISE_FILE_DATE = "13 mai"
FROZEN_RULES_LABEL = "règles avant le 27 mai à 10h-11h"
FROZEN_RULES_REF = "aa78455" # commit ? ou tag ? ou "HEAD" pour la version actuelle (mais moins stable comme référence)

LOCAL_REFERENCE_LABEL = f"fichier Elise du {ELISE_FILE_DATE} + {FROZEN_RULES_LABEL}"


def load_google_sheet(url_or_id: str) -> Tuple[Optional[pd.DataFrame], Optional[str]]:
    if not url_or_id:
        return None, "Empty URL"

    # Extract sheet id
    import re

    match = re.search(r"/d/([A-Za-z0-9_-]+)", url_or_id)
    sheet_id = match.group(1) if match else url_or_id

    # Extract gid if present
    gid_match = re.search(r"gid=([0-9]+)", url_or_id)
    gid = gid_match.group(1) if gid_match else "0"

    csv_url = (
        f"https://docs.google.com/spreadsheets/d/{sheet_id}/export?format=csv&gid={gid}"
    )

    try:
        df = pd.read_csv(csv_url)
        return df, None
    except Exception as exc:
        return None, str(exc)


def load_input(sheet: Optional[str], file_path: Optional[str]) -> pd.DataFrame:
    if sheet and file_path:
        raise ValueError("Use either --sheet or --file, not both.")

    if file_path:
        if file_path.lower().endswith(".csv"):
            return pd.read_csv(file_path)
        if file_path.lower().endswith((".xlsx", ".xls")):
            return pd.read_excel(file_path)
        raise ValueError("Unsupported file type. Use .csv, .xlsx, or .xls.")

    sheet = sheet or DEFAULT_SHEET_URL
    df, err = load_google_sheet(sheet)
    if err:
        raise RuntimeError(f"Failed to load Google Sheet: {err}")
    return df


def input_jobs(sheet: Optional[str], file_path: Optional[str]):
    if sheet and file_path:
        raise ValueError("Use either --sheet or --file, not both.")
    if sheet:
        return [("Google Sheet", lambda: load_input(sheet, None))]
    if file_path:
        return [(os.path.basename(file_path), lambda: load_input(None, file_path))]

    jobs = [("Google Sheet", lambda: load_input(DEFAULT_SHEET_URL, None))]
    if os.path.exists(DEFAULT_LOCAL_FILE):
        jobs.append((DEFAULT_LOCAL_FILE, lambda: load_input(None, DEFAULT_LOCAL_FILE)))
    else:
        print(f"Warning: local file '{DEFAULT_LOCAL_FILE}' not found; skipped.")
    return jobs


def detect_columns(df: pd.DataFrame) -> Tuple[str, Optional[str], Optional[str]]:
    formule_col = None
    for col in df.columns:
        if str(col).strip().lower() == "formule":
            formule_col = col
            break

    if formule_col is None:
        raise ValueError("No 'Formule' column found (case-insensitive).")

    count_i_col = None
    count_j_col = None
    for col in df.columns:
        norm = str(col).strip().lower()
        if norm in ("count_i", "count_iscn"):
            count_i_col = col
        elif norm in ("count_j", "count_jon"):
            count_j_col = col

    return formule_col, count_i_col, count_j_col


def normalize_reference(value):
    if isinstance(value, str):
        trimmed = value.strip()
        if not trimmed:
            return None
        try:
            numeric = float(trimmed)
            if numeric.is_integer():
                return int(numeric)
            return numeric
        except ValueError:
            return trimmed
    if pd.isna(value):
        return None
    if isinstance(value, float) and value.is_integer():
        return int(value)
    return value


def analyze_rows(
    df: pd.DataFrame,
    formule_col: str,
    count_i_col: Optional[str],
    count_j_col: Optional[str],
    analyzer: Callable = analyser_formule,
):
    audit_rows = []
    errors = 0
    match_i_total = 0
    match_i_ok = 0
    match_j_total = 0
    match_j_ok = 0

    for idx, row in df.iterrows():
        formule = row[formule_col]
        ref_i = normalize_reference(row[count_i_col]) if count_i_col else None
        ref_j = normalize_reference(row[count_j_col]) if count_j_col else None

        df_analyse, totals, error = analyzer(formule)
        iscn = None
        jon = None
        if error:
            errors += 1
        else:
            iscn = totals.get("iscn")
            jon = totals.get("jondroville")

        match_i = None
        match_j = None
        if count_i_col and ref_i is not None:
            match_i = iscn == ref_i
            match_i_total += 1
            if match_i:
                match_i_ok += 1
        if count_j_col and ref_j is not None:
            match_j = jon == ref_j
            match_j_total += 1
            if match_j:
                match_j_ok += 1

        audit_rows.append(
            {
                "source": None,
                "index": idx + 1,
                "formule": formule,
                "iscn_obtenu": iscn if error is None else None,
                "jon_obtenu": jon if error is None else None,
                "ref_iscn": ref_i,
                "ref_jon": ref_j,
                "match_iscn": match_i,
                "match_jon": match_j,
                "erreur": error,
            }
        )

    summary = {
        "total": len(df),
        "errors": errors,
        "match_i_total": match_i_total,
        "match_i_ok": match_i_ok,
        "match_j_total": match_j_total,
        "match_j_ok": match_j_ok,
    }
    return audit_rows, summary


def discordant_rows(audit_rows):
    rows = []
    for row in audit_rows:
        if row["erreur"]:
            rows.append(row)
            continue
        if row["match_iscn"] is False or row["match_jon"] is False:
            rows.append(row)
    return rows


def print_summary(summary: dict, label: Optional[str] = None):
    total = summary["total"]
    errors = summary["errors"]
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    if label:
        print(f"📄 Source: {label}")
    print(f"⏱  Execution  à : {now}")
    print(f"    - Total lignes: {total}")
    print(f"    - Erreurs d'analyse: {errors}")

    def color(text: str, code: str) -> str:
        return f"\033[{code}m{text}\033[0m"

    color_iscn = "34"  # blue
    color_jon = "33"   # brown-ish (yellow)

    def format_percent(rate: float) -> str:
        rounded_down = math.floor(rate * 100) / 100
        return f"{rounded_down:.2f}".replace(".", ",")

    def progress_bar(rate: float, width: int = 20) -> str:
        filled = int(math.floor((rate / 100) * width))
        filled = max(0, min(width, filled))
        return f"[{'█' * filled}{'░' * (width - filled)}]"

    if summary["match_i_total"]:
        rate = (summary["match_i_ok"] / summary["match_i_total"]) * 100
        label = f"({color('ISCN', color_iscn)})"
        print(
            f"    - Match ISCN {label}: {summary['match_i_ok']}/{summary['match_i_total']} ({format_percent(rate)}%) {progress_bar(rate)}"
        )
    else:
        print(f"    - Match ISCN ({color('ISCN', color_iscn)}): N/A (aucune reference)")

    if summary["match_j_total"]:
        rate = (summary["match_j_ok"] / summary["match_j_total"]) * 100
        label = f"({color('Jon', color_jon)})"
        print(
            f"    - Match Jon {label}: {summary['match_j_ok']}/{summary['match_j_total']} ({format_percent(rate)}%) {progress_bar(rate)}"
        )
    else:
        print(f"    - Match Jon ({color('Jon', color_jon)}): N/A (aucune reference)")


def format_match(value):
    if value is None:
        return "N/A"
    if value is True:
        return "TRUE"
    return "\033[31mFALSE\033[0m"


def _normalize_index(value):
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except Exception:
        pass
    if isinstance(value, float) and value.is_integer():
        return int(value)
    return value


def _index_map(rows):
    mapped = {}
    for row in rows:
        source = row.get("source")
        idx = _normalize_index(row.get("index"))
        if idx is not None:
            mapped[(source, idx)] = row
    return mapped


def _norm_value(value):
    try:
        if pd.isna(value):
            return None
    except Exception:
        pass
    if isinstance(value, float) and value.is_integer():
        return int(value)
    return value


def _differs(a, b):
    status_a = _status_from_row(a)
    status_b = _status_from_row(b)
    if status_a or status_b:
        if status_a != status_b:
            return True
        return False

    keys = [
        "iscn_obtenu",
        "jon_obtenu",
        "ref_iscn",
        "ref_jon",
        "match_iscn",
        "match_jon",
        "erreur",
    ]
    comparable = False
    for key in keys:
        if key not in a or key not in b:
            continue
        comparable = True
        if _norm_value(a.get(key)) != _norm_value(b.get(key)):
            return True
    if comparable:
        return False
    return False


def print_discordances(
    audit_rows,
    limit=20,
    baseline_rows=None,
    baseline_label=None,
    source_label=None,
):
    discordances = discordant_rows(audit_rows)

    print(f"\n⚠️  Discordances  (nb total = {len(discordances)}) :")

    def color(text: str, code: str) -> str:
        return f"\033[{code}m{text}\033[0m"

    color_iscn = "34"
    color_jon = "33"

    def format_index(row):
        if row["erreur"]:
            return f"{row['index']}(Erreur)"
        iscn_bad = row["match_iscn"] is False
        jon_bad = row["match_jon"] is False
        if iscn_bad and jon_bad:
            return str(row["index"])
        if iscn_bad:
            return f"{row['index']}({color('ISCN', color_iscn)})"
        if jon_bad:
            return f"{row['index']}({color('Jon', color_jon)})"
        return str(row["index"])

    shown_indices = [format_index(row) for row in discordances]
    print(f"Indices: {', '.join(shown_indices) if shown_indices else '—'}")

    if baseline_rows is not None:
        baseline_rows = _baseline_rows_for_source(baseline_rows, source_label)
        if baseline_label:
            print(f"\nComparaison baseline: {baseline_label}")
        base_map = _index_map(baseline_rows)
        new_map = _index_map(discordances)

        new_indices = [format_index(new_map[key]) for key in new_map.keys() if key not in base_map]
        resolved = [str(key[1]) for key in base_map.keys() if key not in new_map]
        perturbed = [
            str(key[1])
            for key in new_map.keys()
            if key in base_map and _differs(new_map[key], base_map[key])
        ]

        print(f"Nouvelles: {', '.join(new_indices) if new_indices else '—'}")
        print(f"Corrigées: {', '.join(resolved) if resolved else '—'}")
        print(f"Différentes: {', '.join(perturbed) if perturbed else '—'}")

    return discordances


def _baseline_rows_for_source(rows, source_label):
    if not source_label:
        return rows

    has_source = any(row.get("source") for row in rows)
    if not has_source:
        return rows
    return [row for row in rows if row.get("source") == source_label]


def _load_baseline_csv(path: str):
    try:
        baseline_df = pd.read_csv(path)
        return baseline_df.to_dict(orient="records")
    except pd.errors.EmptyDataError:
        return None
    except Exception as exc:
        print(f"Warning: unable to load baseline '{path}': {exc}")
        return None


def _load_baseline_from_git(path: str):
    try:
        git_root = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"],
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
        rel_path = os.path.relpath(os.path.abspath(path), git_root)
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
        content = subprocess.check_output(
            ["git", "show", f"HEAD:{rel_path}"],
            stderr=subprocess.DEVNULL,
            text=True,
        )
        baseline_df = pd.read_csv(io.StringIO(content))
        return baseline_df.to_dict(orient="records"), commit
    except Exception:
        return None, None


def _load_analyzer_from_git(ref: str = "HEAD"):
    try:
        source = subprocess.check_output(
            ["git", "show", f"{ref}:My_expert_karyo_functions.py"],
            stderr=subprocess.DEVNULL,
            text=True,
        )
        module = types.ModuleType(f"baseline_my_expert_{ref.replace('/', '_')}")
        exec(compile(source, f"{ref}:My_expert_karyo_functions.py", "exec"), module.__dict__)
        return module.analyser_formule
    except Exception as exc:
        print(f"Warning: unable to load baseline analyzer from {ref}: {exc}")
        return None


def _baseline_rows_from_git_analyzer(
    df: pd.DataFrame,
    formule_col: str,
    count_i_col: Optional[str],
    count_j_col: Optional[str],
    source_label: str,
):
    analyzer = _load_analyzer_from_git(FROZEN_RULES_REF)
    if analyzer is None:
        return None

    rows, _ = analyze_rows(df, formule_col, count_i_col, count_j_col, analyzer=analyzer)
    for row in rows:
        row["source"] = source_label
    return discordant_rows(rows)


def _is_default_local_source(label: str) -> bool:
    return os.path.basename(label) == DEFAULT_LOCAL_FILE


def _format_status(label: str, obtenu, ref, match):
    if ref is None:
        return f"{label} N/A"
    if match is True:
        return f"{label} ok {obtenu}/{ref}"
    if match is False:
        return f"{label} False {obtenu} au lieu de {ref}"
    return f"{label} N/A"


def _status_from_row(row):
    def _clean(value):
        try:
            if pd.isna(value):
                return ""
        except Exception:
            pass
        return str(value)

    if "iscn_status" in row or "jon_status" in row:
        return (
            _clean(row.get("iscn_status", "")),
            _clean(row.get("jon_status", "")),
            _clean(row.get("erreur", "")),
        )
    return (
        _format_status("ISCN", row.get("iscn_obtenu"), row.get("ref_iscn"), row.get("match_iscn")),
        _format_status("Jon", row.get("jon_obtenu"), row.get("ref_jon"), row.get("match_jon")),
        _clean(row.get("erreur", "")),
    )


def export_audit_csv(audit_rows, out_path: str):
    rows = []
    for row in audit_rows:
        is_discordant = bool(
            row.get("erreur")
            or row.get("match_iscn") is False
            or row.get("match_jon") is False
        )
        rows.append(
            {
                "source": row.get("source"),
                "index": row.get("index"),
                "formule": row.get("formule"),
                "iscn_status": _format_status(
                    "ISCN",
                    row.get("iscn_obtenu"),
                    row.get("ref_iscn"),
                    row.get("match_iscn"),
                ),
                "jon_status": _format_status(
                    "Jon",
                    row.get("jon_obtenu"),
                    row.get("ref_jon"),
                    row.get("match_jon"),
                ),
                "erreur": row.get("erreur") or "",
            }
        )
    columns = ["source", "index", "formule", "iscn_status", "jon_status", "erreur"]
    df_audit = pd.DataFrame(rows, columns=columns)
    df_audit.to_csv(out_path, index=False)


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Compare analyser_formule results against a sheet or local file."
    )
    parser.add_argument(
        "--sheet",
        help="Google Sheet URL or id (defaults to built-in sheet if omitted).",
    )
    parser.add_argument(
        "--file",
        help="Local CSV/XLSX file path (exclusive with --sheet).",
    )
    parser.add_argument(
        "--out",
        default="audit.csv",
        help="Output CSV audit path (default: audit.csv).",
    )
    parser.add_argument(
        "--baseline",
        help="Previous audit CSV to flag new discordances.",
    )
    parser.add_argument(
        "--prev",
        action="store_true",
        help="Compare only against previous audit CSV instead of commit baseline.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=20,
        help="Number of discordances to display (default: 20).",
    )
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv or sys.argv[1:])
    try:
        baseline_rows = None
        baseline_label = None
        baseline_path = args.baseline
        if baseline_path is None and os.path.exists(args.out):
            baseline_path = args.out
        if baseline_path:
            baseline_rows = _load_baseline_csv(baseline_path)
            if baseline_rows is not None:
                baseline_label = os.path.basename(baseline_path)

        commit_rows, commit_ref = _load_baseline_from_git(args.out)
        all_discordances = []
        jobs = input_jobs(args.sheet, args.file)
        for job_index, (label, loader) in enumerate(jobs):
            if job_index:
                print("\n" + "=" * 72 + "\n")

            df = loader()
            formule_col, count_i_col, count_j_col = detect_columns(df)
            audit_rows, summary = analyze_rows(df, formule_col, count_i_col, count_j_col)
            for row in audit_rows:
                row["source"] = label

            run_baseline_rows = baseline_rows
            run_baseline_label = baseline_label
            if (
                _is_default_local_source(label)
                and args.baseline is None
                and not args.prev
            ):
                old_rows = _baseline_rows_from_git_analyzer(
                    df,
                    formule_col,
                    count_i_col,
                    count_j_col,
                    label,
                )
                if old_rows is not None:
                    run_baseline_rows = old_rows
                    run_baseline_label = LOCAL_REFERENCE_LABEL

            print_summary(summary, label=label)
            if args.prev or run_baseline_rows is not None:
                discordances = print_discordances(
                    audit_rows,
                    limit=args.limit,
                    baseline_rows=run_baseline_rows,
                    baseline_label=run_baseline_label,
                    source_label=label,
                )
            else:
                discordances = print_discordances(
                    audit_rows,
                    limit=args.limit,
                    baseline_rows=commit_rows,
                    baseline_label=f"commit {commit_ref[:8]}" if commit_ref else None,
                    source_label=label,
                )
            all_discordances.extend(discordances)

        if os.path.exists(args.out):
            prev_path = f"{args.out}.prev"
            try:
                shutil.copyfile(args.out, prev_path)
            except Exception as exc:
                print(f"Warning: unable to save previous audit to '{prev_path}': {exc}")
        export_audit_csv(all_discordances, args.out)
        print(f"Audit CSV ecrit: {args.out}")
        return 0
    except Exception as exc:
        print(f"Erreur: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
