#!/usr/bin/env python3
"""Régénère le CSV interne du référentiel à partir du catalogue Python."""

import ast
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "My_expert_karyo_functions.py"
DESTINATION = ROOT / "rules_catalog_reference.csv"

# Textes anglais canoniques. Les mêmes noms de colonnes doivent être utilisés
# dans le Google Sheet public pour permettre la comparaison automatique.
ENGLISH_TEXTS = {
    "ISCN.ABSENT_IN_REPEAT_ZERO": ("Absent in a repeated clone", "An anomaly explicitly absent from a secondary clone and already counted previously."),
    "ISCN.ABSENT_IN_REPEAT_IMPLICIT_LOSS": ("Implicit loss in a repeated clone", "Absence of a structural anomaly interpreted as an implicit chromosome loss."),
    "ISCN.CONSTITUTIONAL_GAIN": ("Constitutional gain", "Chromosome gain annotated as constitutional."),
    "ISCN.CONSTITUTIONAL_SUSPECT": ("Suspected constitutional anomaly", "Anomaly annotated as possibly constitutional."),
    "ISCN.CONSTITUTIONAL_CERTAIN": ("Confirmed constitutional anomaly", "Anomaly annotated as constitutional."),
    "ISCN.REPEAT_NOTATION": ("Repeat notation", "idem/sl/sdl notation already represented by a previous clone."),
    "ISCN.INTERTWINED_DER_BALANCED_SIMPLE": ("Balanced intertwined derivatives", "Group of related derivatives with concordant breakpoints."),
    "ISCN.INTERTWINED_DER_BALANCED_INSERTION": ("Balanced intertwined derivatives with insertion", "Balanced group of related derivatives with an additional insertion."),
    "ISCN.INTERTWINED_DER_UNBALANCED": ("Unbalanced intertwined derivatives", "Group of related derivatives without sufficient breakpoint concordance for a balanced event."),
    "ISCN.INTERTWINED_DER_PART": ("Derivative already included in an intertwined event", "Derivative belonging to an intertwined rearrangement already counted."),
    "ISCN.BALANCED_T": ("Balanced translocation", "Translocation with compatible reciprocal derivatives."),
    "ISCN.UNBALANCED_T": ("Unbalanced translocation", "Explicitly detected unbalanced translocation."),
    "ISCN.T_ALREADY_IN_DER": ("Translocation already counted", "Explicit translocation already included in a derivative chromosome."),
    "ISCN.IMPLICIT": ("Implicit anomaly", "Anomaly already explained by a reference anomaly."),
    "ISCN.STRUCTURAL_GAIN_DUPLICATE": ("Gain of an existing structural anomaly", "Additional copy of a structural anomaly already described."),
    "ISCN.SEMANTIC_PLUS_ISO": ("Chromosome gain with isochromosome", "The +i(...) notation combines a chromosome gain and an isochromosome."),
    "ISCN.SEMANTIC_PLUS_DEL": ("Chromosome gain with deletion", "The +del(...) notation combines a chromosome gain and a deletion."),
    "ISCN.MAR": ("Marker chromosome", "Marker chromosome, including quantified forms."),
    "ISCN.DICENTRIC": ("Dicentric chromosome", "Anomaly identified as a dicentric chromosome."),
    "ISCN.JUMPING_TRANSLOCATION_SECONDARY": ("Secondary jumping translocation occurrence", "Secondary derivative sharing the same donor chromosome and breakpoint in another subclone."),
    "ISCN.T_VIA_DER_BALANCED": ("Balanced translocation via derivative", "Balanced translocation counted from a derivative chromosome."),
    "ISCN.T_VIA_DER_UNBALANCED": ("Unbalanced translocation via derivative", "Unbalanced translocation counted from a derivative chromosome."),
    "ISCN.DER_T_COUNTED": ("Associated derivative already counted", "Derivative associated with a translocation that has already been counted."),
    "ISCN.DER_BALANCED_T": ("Balanced-translocation derivative", "Derivative originating from an already represented balanced translocation."),
    "ISCN.DER_NO_BREAKPOINT": ("Derivative without a detailed breakpoint", "Derivative chromosome without a detailed breakpoint."),
    "ISCN.DER_MULTI": ("Multi-chromosome derivative", "Derivative involving several identified chromosomes."),
    "ISCN.DER_MULTI_UNCERTAIN": ("Uncertain multi-chromosome derivative", "Derivative involving several chromosomes with uncertain notation."),
    "ISCN.DER_UNCERTAIN_SECOND": ("Derivative with an uncertain second chromosome", "Derivative that may involve a second chromosome that cannot be identified with certainty."),
    "ISCN.DER_SAME_CHR": ("Same-chromosome derivative", "Derivative involving a single chromosome or an intrachromosomal rearrangement."),
    "ISCN.SINGLE_CHR_GAIN_REPEAT": ("Repeated single-chromosome gain", "Repeated gain of the same chromosome."),
    "ISCN.SINGLE_CHR_TRIPLICATION": ("Triplication", "Triplication involving a single chromosome."),
    "ISCN.SINGLE_CHR_ISODERIVATIVE": ("Isoderivative chromosome", "Isoderivative or isodicentric chromosome."),
    "ISCN.COMPLEX_MULTI_CHR": ("Complex multi-chromosome anomaly", "Unbalanced anomaly involving several chromosomes."),
    "ISCN.UNBALANCED_TRANSLOCATION": ("Unbalanced translocation", "Non-pure translocation or translocation carried by a derivative."),
    "ISCN.GAIN_SIMPLE": ("Simple gain", "Simple chromosome gain."),
    "ISCN.LOSS_SIMPLE": ("Simple loss", "Simple chromosome loss."),
    "ISCN.OTHER_STANDARD": ("Other standard anomaly", "Non-constitutional anomaly not covered by a more specific rule."),
    "JON.ABSENT_IN_REPEAT_ZERO": ("Absent in a repeated clone", "An anomaly explicitly absent from a secondary clone and already counted previously."),
    "JON.ABSENT_IN_REPEAT_IMPLICIT_LOSS": ("Implicit loss in a repeated clone", "Absence of a structural anomaly interpreted as an implicit chromosome loss."),
    "JON.CONSTITUTIONAL_GAIN": ("Constitutional gain", "Chromosome gain annotated as constitutional."),
    "JON.CONSTITUTIONAL_SUSPECT": ("Suspected constitutional anomaly", "Anomaly annotated as possibly constitutional."),
    "JON.CONSTITUTIONAL_CERTAIN": ("Confirmed constitutional anomaly", "Anomaly annotated as constitutional."),
    "JON.REPEAT_NOTATION": ("Repeat notation", "idem/sl/sdl notation already represented by a previous clone."),
    "JON.IMPLICIT": ("Implicit anomaly", "Duplication of a reference anomaly."),
    "JON.MAR": ("Marker chromosome", "Marker chromosome."),
    "JON.TRIPLOIDY": ("Triploidy", "Triploidy excluded from the Jondreville count."),
    "JON.DEFAULT": ("Non-constitutional anomaly", "Each non-constitutional anomaly counts as one point."),
}


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
            "Libellé v_EN": ENGLISH_TEXTS[rule_id][0],
            "Explication v_EN": ENGLISH_TEXTS[rule_id][1],
            "Critère technique v_EN": (
                f"Detection criterion: {ENGLISH_TEXTS[rule_id][1]} "
                f"Default score: {score}."
            ),
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
        "Libellé v_EN",
        "Explication v_EN",
        "Critère technique v_EN",
    ]

    with DESTINATION.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"{DESTINATION.name}: {len(rows)} règles écrites")


if __name__ == "__main__":
    main()
