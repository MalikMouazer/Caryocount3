#!/usr/bin/env python3
"""Régénère le CSV interne du référentiel à partir du catalogue Python."""

import ast
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "My_expert_karyo_functions.py"
DESTINATION = ROOT / "rules_catalog_reference.csv"


def assigned_value(tree: ast.Module, name: str) -> ast.AST:
    for node in tree.body:
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            if any(isinstance(target, ast.Name) and target.id == name for target in targets):
                return node.value
    raise RuntimeError(f"Affectation {name} introuvable dans {SOURCE.name}")


def main() -> None:
    tree = ast.parse(SOURCE.read_text(encoding="utf-8"), filename=str(SOURCE))
    catalog_node = assigned_value(tree, "RULE_CATALOG")
    checks = ast.literal_eval(assigned_value(tree, "RULE_TECHNICAL_CHECKS"))
    priorities = ast.literal_eval(assigned_value(tree, "RULE_PRIORITY"))

    rows_by_id = {}
    for item in catalog_node.elts:
        if not isinstance(item, ast.Call) or not isinstance(item.func, ast.Name):
            continue
        if item.func.id != "RuleSpec":
            continue
        rule_id, system, score, title, explanation = (
            ast.literal_eval(argument) for argument in item.args
        )
        rows_by_id[rule_id] = {
            "Rule ID": rule_id,
            "Référentiel": system,
            "Score par défaut": score,
            "Critère technique": checks.get(rule_id, ""),
            "Libellé": title,
            "Explication": explanation,
        }

    ordered_ids = [
        rule_id
        for system in ("ISCN 2024", "Jondreville 2020")
        for rule_id in priorities.get(system, ())
    ]
    rows = [rows_by_id[rule_id] for rule_id in ordered_ids]
    counters = {"ISCN 2024": 0, "Jondreville 2020": 0}
    for row in rows:
        system = row["Référentiel"]
        counters[system] += 1
        suffix = "JON" if system == "Jondreville 2020" else "ISCN"
        row["N°"] = f"{counters[system]}_{suffix}"

    fieldnames = [
        "Rule ID",
        "N°",
        "Référentiel",
        "Score par défaut",
        "Critère technique",
        "Libellé",
        "Explication",
    ]

    with DESTINATION.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"{DESTINATION.name}: {len(rows)} règles écrites")


if __name__ == "__main__":
    main()
