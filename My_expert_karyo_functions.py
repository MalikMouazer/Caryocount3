import re
from collections import Counter
import pandas as pd
from dataclasses import dataclass
import os
from urllib.parse import parse_qs, urlparse

# =========================
# Parsing
# =========================

# Extraction des numéros de chromosome dans une anomalie ISCN
def get_chromosomes(anom):
    """Retourne l'ensemble des chromosomes impliqués dans ``anom``.

    La fonction détecte les chromosomes apparaissant :
    - juste après les mots clés (der, del, dup, t, ...)
    - dans la seconde parenthèse des notations ``der(...)`` (apès les flèches)
    - précédés d'un ``?`` comme dans ``t(?1;17)``
    """

    nums: set[str] = set()

    # 1) Chromosomes directement après der(...), t(...), etc.
    for m in re.finditer(r'(?:der|dic|del|dup|ins|t|i|ider|idic|r)\((\??[0-9XYxy;?]+)', anom):
        raw = m.group(1)
        # Se limiter à la partie chromosomique avant un ")" ou une nouvelle parenthèse
        raw = re.split(r'[)()]', raw)[0]
        for num in raw.split(';'):
            cleaned = num.lstrip('?').upper()
            if cleaned:
                nums.add(cleaned)
            elif '?' in num:
                nums.add('?')

    # 2) Chromosomes mentionnés dans la seconde parenthèse des der(...)
    for _, second in re.findall(r'der\(([^)]*)\)\(([^)]*)\)', anom):
        for n in re.findall(r'\??(\d+|X|Y)(?=[pq])', second, re.IGNORECASE):
            nums.add(n.lstrip('?').upper())
        if '?' in second:
            nums.add('?')

    return nums


def count_known_chromosomes(chroms: set[str]) -> int:
    """Compte uniquement les chromosomes identifiés (sans '?')."""

    return sum(1 for c in chroms if c != '?')

# Parsing de la formule karyotypique
def expand_condensed_clone(clone: str) -> tuple[str, bool]:
    """Déplie une écriture condensée de type 50(25,X,+X,+21)×2[9]."""

    pattern = re.compile(r"^(\d+)\(([^)]*)\)[x×](\d+)(\[[^\]]+\])?$")
    m = pattern.match(clone)
    if not m:
        return clone, False

    total, inner, mult_raw, suffix = m.groups()
    mult = int(mult_raw)
    tokens = [t for t in inner.split(",") if t]
    sex_x = 0
    sex_y = 0
    expanded: list[str] = []
    for tok in tokens:
        if re.fullmatch(r"\d+", tok):
            continue
        if tok in ("X", "Y"):
            if tok == "X":
                sex_x += mult
            else:
                sex_y += mult
            continue
        expanded.extend([tok] * mult)

    sex_token = ""
    if sex_x or sex_y:
        sex_token = ("X" * sex_x) + ("Y" * sex_y)

    parts = []
    if sex_token:
        parts.append(sex_token)
    parts.extend(expanded)

    expanded_str = ",".join(parts)
    base = f"{total}" if not expanded_str else f"{total},{expanded_str}"
    return f"{base}{suffix or ''}", True


def expand_condensed_formula(formule: str) -> tuple[str, bool]:
    """Déplie les notations condensées dans une formule ISCN."""

    cleaned = re.sub(r"\s+", "", formule)
    clones = cleaned.split("/")
    changed = False
    expanded_clones = []
    for clone in clones:
        expanded, did_change = expand_condensed_clone(clone)
        changed = changed or did_change
        expanded_clones.append(expanded)
    return "/".join(expanded_clones), changed


def parse_caryotype(chaine_iscn):
    """Parse une chaîne ISCN et renvoie les anomalies par clone.

    Returns
    -------
    tuple[list[str], dict[str, list[str]], list[dict[str, str]]]
        - Liste plate des anomalies (ordre d'apparition)
        - Mapping ``anomalie -> [clone1, clone2, ...]``
        - Liste d'entrées ``{"anomaly": str, "clone": str}``
    """
    # Remove all whitespace for robust parsing
    chaine_iscn = re.sub(r"\s+", "", chaine_iscn)

    anomalies = []
    clone_map = {}
    entries: list[dict[str, str]] = []
    clones_raw = chaine_iscn.split('/')
    clones = [expand_condensed_clone(c)[0] for c in clones_raw]
    clones = [re.sub(r"\[.*?\]", "", c) for c in clones]
    for idx, clone in enumerate(clones, start=1):
        clone_name = f"clone{idx}"
        parts = [p.strip().strip('.') for p in clone.split(',') if p.strip()]
        clone_has_repeat = any(is_repeat_notation(p) for p in parts)
        # Détection de la ploidie
        try:
            ploidy_token = parts[0]
            # Les notations XX<2n> servent à indiquer les métaphases et ne doivent
            # pas être confondues avec une anomalie de ploidie.
            if '<2n>' not in ploidy_token:
                match = re.search(r"\d+", ploidy_token)
                total = int(match.group()) if match else None
            else:
                total = None

            if total and total != 46:
                if 58 <= total <= 80:
                    pl = 'Triploidy'
                elif 81 <= total <= 103:
                    pl = 'Tetraploidy'
                else:
                    pl = None

                if pl:
                    anomalies.append(pl)
                    clone_map.setdefault(pl, []).append(clone_name)
                    entries.append({"anomaly": pl, "clone": clone_name})
        except Exception:
            pass
        # Extraction des anomalies structurelles
        for an in parts[2:]:
            count = anomaly_occurrences(an)
            anomalies.extend([an] * count)
            clone_map.setdefault(an, []).append(clone_name)
            entries.append({
                "anomaly": an,
                "clone": clone_name,
                "count": count,
                "repeat_clone": clone_has_repeat,
            })
    return anomalies, clone_map, entries


def display_clone_label(clone_name: str) -> str:
    """Retourne un libellé "Clone N" lisible pour un identifiant interne."""

    if not clone_name:
        return ""
    m = re.match(r"clone(\d+)", clone_name, re.IGNORECASE)
    if m:
        return f"Clone {m.group(1)}"
    return clone_name.capitalize()


# =========================
# Typage / Détection
# =========================
def is_repeat_notation(anom: str) -> bool:
    """Détecte les notations de répétition d'anomalies (idem, sl, sdl...)."""
    base = anom.strip().lower()
    if base == 'idem':
        return True
    if re.match(r'^(?:sl|sdl)\d*$', base):
        return True
    return False

# Détection des anomalies unichromosomiques déséquilibrées de poids 2
def is_single_chr_deseq(anom, count):
    """
    Détecte les anomalies unichromosomiques déséquilibrées qui valent 2 points:
    - Tetrasomie/triplication/quadruplication
    - Chromosome isodérivé
    """
    # Tetrasomie/triplication/quadruplication
    if anom.startswith('+') and count > 1:
        return True
    if anom.startswith('trp'):
        return True
    # Chromosome isodérivé ou isodicentrique
    if anom.startswith('ider'):
        return True
    return False

# Détection des anomalies équilibrées
def is_balanced_translocation(anom):
    """
    Détecte les translocations équilibrées:
    t(NUM;NUM[;...])(p;q) sans der,+,-
    """
    pattern = r'^t\(\??(?:\d+|[XY])(?:;\??(?:\d+|[XY]))+\)\(.+\)$'
    if 'der' in anom or '+' in anom:
        return False

    match = re.match(pattern, anom, re.IGNORECASE)
    if not match:
        return False

    # Retirer la partie décrivant les points de cassure pour ignorer
    # les tirets utilisés dans les plages (ex: q12-13)
    breakpoint_start = anom.rfind('(')
    non_breakpoint_part = anom[:breakpoint_start] if breakpoint_start != -1 else anom

    if '-' in non_breakpoint_part:
        return False

    return True

def is_unbalanced_translocation(anom):
    """
    Détecte les translocations déséquilibrées:
    - chromosome dérivé (der(...)) contenant un t(...) ou
    - tout t(...) non pure
    """
    # Cas d'un chromosome dérivé ou dicentrique comportant une translocation
    if ('der' in anom or 'dic' in anom) and 't(' in anom:
        return True
    # Cas d'un t(...) quelconque non pur (équilibré)
    if 't(' in anom and not is_balanced_translocation(anom):
        return True
    return False

def is_balanced_insertion(anom):
    """
    Détecte les insertions équilibrées:
    ins(NUM;NUM[;...])(p;q1q2) sans der,+,-
    """
    pattern = r'^ins\(\??\d+(?:;\??\d+)+\)\(.+\)$'
    return bool(re.match(pattern, anom)) and 'der' not in anom and '+' not in anom and '-' not in anom

# Détection des anomalies multichromosomiques déséquilibrées pour 2 points
def is_complex_multichr_deseq(anom):
    """
    Détecte les anomalies multichromosomiques déséquilibrées (≥2 chromosomes) pour 2 points.
    Les chromosomes dérivés ("der") sans points de cassure détaillés ne sont
    pas considérés complexes: on ne sait pas s'ils impliquent un ou deux
    chromosomes.
    Renvoie False si un seul chromosome impliqué.
    """
    # Cas particulier des chromosomes dérivés
    if anom.startswith('der'):
        if is_der_without_breakpoints(anom):
            return False
        chroms = get_chromosomes(anom)
        known = count_known_chromosomes(chroms)
        # s'il n'y a pas au moins deux chromosomes identifiés -> pas complexe
        return known >= 2

    chroms = get_chromosomes(anom)
    known = count_known_chromosomes(chroms)
    # si un seul chromosome identifié -> pas multi-chromosomique déséquilibrée
    if known <= 1:
        return False
    # chromosome dicentrique ou anneau -> complexe multi-chromosomique
    if anom.startswith('dic') or anom.startswith('r('):
        return True
    # insertion non pure -> complexe
    if 'ins(' in anom and not is_balanced_insertion(anom):
        return True
    # translocation non pure -> complexe
    if 't(' in anom and not is_balanced_translocation(anom):
        return True
    return False

# Typage pour affichage
def type_anomalie(anom):
    """
    Détermine le type d'anomalie pour l'affichage.
    Retourne une chaîne décrivant le type d'anomalie.
    """
    base = strip_multiplicity(strip_sign(anom))
    # Si une anomalie structurale est précédée d'un '+', ne pas l'interpréter
    # comme un gain simple de chromosome.
    if anom.startswith("+"):
        if base.startswith("del"):
            return "Délétion"
        if base.startswith("dup"):
            return "Duplication"
        if base.startswith("inv"):
            return "Inversion"
        if base.startswith("ins"):
            return "Insertion"
        if base.startswith("add"):
            return "Addition"
        if base.startswith("der"):
            return "Chromosome dérivé"
        if base.startswith("dic"):
            return "Chromosome dicentrique"
        if base.startswith("idic"):
            return "Isodicentric chromosome"
        if base.startswith("ider"):
            return "Isoderivative chromosome"
        if base.startswith("i(") or "iso" in base:
            return "Isochromosome"
        if base.startswith("t("):
            return "Translocation"
    if is_negative_anomaly(anom):
        if base.startswith("del"):
            return "Délétion"
        if base.startswith("dup"):
            return "Duplication"
        if base.startswith("inv"):
            return "Inversion"
        if base.startswith("ins"):
            return "Insertion"
        if base.startswith("add"):
            return "Addition"
        if base.startswith("der"):
            return "Chromosome dérivé"
        if base.startswith("dic"):
            return "Chromosome dicentrique"
        if base.startswith("idic"):
            return "Isodicentric chromosome"
        if base.startswith("ider"):
            return "Isoderivative chromosome"
        if base.startswith("i(") or "iso" in base:
            return "Isochromosome"
        if base.startswith("t("):
            return "Translocation"
    if is_marker_anomaly(anom):
        m = re.match(r'^mar(\d+)$', base)
        suffix = m.group(1) if m else ""
        return f"Gain mar{suffix}"
    if is_repeat_notation(anom):
        return 'Notation de répétition d’anomalies'
    if anom == '<2n>':
        return 'Ploidy'
    if '~' in anom:
        return 'Pléiade chromosomique'
    base = strip_sign(anom)
    if base.startswith('mar'):
        return f"Gain {base}"
    if 'dmin' in anom:
        return 'Double minutes'
    if anom.startswith('hsr'):
        return 'Homogeneously staining region'
    if anom.startswith('r('):
        return 'Anneau'
    if anom.startswith('der'):
        return 'Chromosome dérivé'
    if anom.startswith('ins'):
        return 'Insertion'
    if anom.startswith('inv'):
        return 'Inversion'
    if anom.startswith('t('):
        return 'Translocation'
    if anom.startswith('add'):
        return 'Addition'
    if anom == 'Triploidy':
        return 'Triploïdie'
    if anom == 'Tetraploidy':
        return 'Tétraploïdie'
    if anom.startswith('+'):
        return 'Gain chr' + re.sub(r"\D", "", anom)
    if anom.startswith('-'):
        return 'Perte chr' + re.sub(r"\D", "", anom)
    if anom.startswith('dup'):
        return 'Duplication'
    if anom.startswith('del'):
        return 'Délétion'
    if anom.startswith('trp'):
        return 'Triplication/Quadruplication'
    if anom.startswith('dic'):
        return 'Chromosome dicentrique'
    if anom.startswith('idic'):
        return 'Isodicentric chromosome'
    if anom.startswith('ider'):
        return 'Isoderivative chromosome'
    if anom.startswith('i(') or 'iso' in anom:
        return 'Isochromosome'
    if is_complex_multichr_deseq(anom):
        return 'Multichromosomique déséquilibrée'
    if is_balanced_translocation(anom):
        return 'Translocation équilibrée'
    if is_unbalanced_translocation(anom):
        return 'Translocation déséquilibrée'
    if is_balanced_insertion(anom):
        return 'Insertion équilibrée'
    return '-'

# Calcul des scores
# =========================
# Normalisation
# =========================
def normalize_anomaly(anom: str) -> str:
    """Normalise une anomalie pour le scoring.

    - Supprime un éventuel point d'interrogation en début d'anomalie
      ("?dic" -> "dic").
    """
    norm = anom.lstrip('?')
    return norm


def strip_sign(anom: str) -> str:
    """Retire un signe initial + ou - sans modifier le reste."""

    return anom[1:] if anom.startswith(("+", "-", "−")) else anom


def is_negative_anomaly(anom: str) -> bool:
    """Indique si l'anomalie commence par un signe négatif (ASCII ou Unicode)."""

    return anom.startswith(("-", "−"))


def strip_negative(anom: str) -> str:
    """Retire un signe négatif initial (ASCII ou Unicode) sans modifier le reste."""

    return anom[1:] if is_negative_anomaly(anom) else anom


def strip_multiplicity(anom: str) -> str:
    """Retire un suffixe de multiplicité (xN ou ×N) sans modifier le reste."""

    return re.sub(r"(?:x|×)\d+$", "", anom)


def is_der_without_breakpoints(anom: str) -> bool:
    """Indique si ``anom`` est un der(x) sans points de cassure détaillés."""

    return bool(re.fullmatch(r"der\([^)]+\)", anom, re.IGNORECASE))


def is_marker_anomaly(anom: str) -> bool:
    """Indique si l'anomalie correspond à un marqueur (mar, +mar, +Nmar, +A~Bmar)."""

    base = strip_multiplicity(strip_sign(anom)).lower()
    return base.startswith("mar") or base.endswith("mar")


def marker_occurrences(anom: str) -> int:
    """Nombre d'occurrences implicites pour les marqueurs."""

    raw = anom.strip()
    if not raw.startswith("+"):
        return 1
    text = raw.replace("∼", "~")
    low = text.lower()

    # +Nmar
    m = re.match(r"^\+(\d+)mar$", low)
    if m:
        return int(m.group(1))

    # +A~Bmar -> B
    m = re.match(r"^\+(\d+)~(\d+)mar$", low)
    if m:
        return int(m.group(2))

    # +marN or +marNxM (N = identifiant)
    m = re.match(r"^\+mar(\d+)(?:[x×](\d+))?$", low)
    if m:
        mult = int(m.group(2)) if m.group(2) else 1
        return mult

    return 1


def anomaly_occurrences(anom: str) -> int:
    """Nombre d'occurrences implicites pour une anomalie (×N ou marqueurs)."""

    if is_marker_anomaly(anom):
        return marker_occurrences(anom)
    m = re.search(r"(?:x|×)(\d+)$", anom.strip(), re.IGNORECASE)
    if m:
        return int(m.group(1))
    return 1


def effective_occurrences(anom: str, count: int, clone_map: dict[str, list[str]]) -> int:
    """Nombre d'occurrences à scorer.

    Les suffixes xN/×N sont toujours traités comme une multiplicité locale,
    indépendamment du contexte de ploidie.
    """

    return count


def constitutional_status(norm: str) -> tuple[bool, str]:
    """Indique si ``norm`` correspond à une anomalie constitutionnelle."""

    cleaned = norm.strip().lower()
    if re.match(r"^\+\d+c$", cleaned):
        return True, "Gain constitutionnel"
    if cleaned.endswith("?c"):
        return True, "Suspicion d'anomalie constitutionnelle (?c)"
    if cleaned.endswith('c'):
        return True, "Anomalie constitutionnelle"
    return False, ""


def constitutional_rule_decision(system_prefix: str, norm: str) -> "RuleDecision | None":
    """Retourne la règle constitutionnelle précise si ``norm`` en contient une."""

    cleaned = norm.strip().lower()
    if re.match(r"^\+\d+c$", cleaned):
        return RuleDecision(
            rule_id=f"{system_prefix}.CONSTITUTIONAL_GAIN",
            score=0,
            explanation="Gain constitutionnel",
        )
    if cleaned.endswith("?c"):
        return RuleDecision(
            rule_id=f"{system_prefix}.CONSTITUTIONAL_SUSPECT",
            score=0,
            explanation="Suspicion d'anomalie constitutionnelle (?c)",
        )
    if cleaned.endswith("c"):
        return RuleDecision(
            rule_id=f"{system_prefix}.CONSTITUTIONAL_CERTAIN",
            score=0,
            explanation="Anomalie constitutionnelle",
        )
    return None


def append_uncertainty_note(anom: str, explanation: str) -> str:
    """Ajoute une mention d'imprécision si la notation contient un '?'."""

    if '?' not in anom:
        return explanation

    note = "Positions incertaines ('?')"
    if explanation:
        return f"{explanation} — {note}"
    return note


# =========================
# Anomalies implicites
# =========================
def detect_implicit_anomalies(anomalies, clone_map=None):
    """Détecte les anomalies implicites et renvoie un dict.

    Le dict a pour clé l'anomalie normalisée et pour valeur un
    dictionnaire avec la clef ``reason`` décrivant la cause et ``ref``
    l'anomalie de référence à afficher entre parenthèses.
    """
    unique_anoms = list(dict.fromkeys(anomalies))
    norm_counts = Counter(normalize_anomaly(a) for a in unique_anoms)
    # mappage normalisé -> version originale pour l'affichage
    norm_to_orig: dict[str, str] = {}
    first_index: dict[str, int] = {}
    for idx, a in enumerate(unique_anoms):
        norm = normalize_anomaly(a)
        norm_to_orig.setdefault(norm, a)
        first_index.setdefault(norm, idx)

    norm_to_clones: dict[str, set[str]] = {}
    if clone_map:
        for raw_anom, clones in clone_map.items():
            norm_to_clones.setdefault(normalize_anomaly(raw_anom), set()).update(clones)

    implicit = {}

    # 1) Dérivés implicites s'il existe une version explicite (add/del/dup)
    t_events = {}
    for an in norm_counts:
        m = re.match(r"(?:der|dic)\((\d+)\).*t\((\d+);(\d+)\)", an)
        if m:
            _, A, B = m.groups()
            key = tuple(sorted([A, B]))
            t_events.setdefault(key, []).append(an)
    for ders in t_events.values():
        explicits = [d for d in ders if re.search(r"add|del|dup", d)]
        if explicits:
            ref = norm_to_orig[explicits[0]]
            for d in ders:
                if d not in explicits:
                    implicit[d] = {"reason": "Dérivé implicite", "ref": ref}

    # 2) Gains/pertes simples issus d'un dérivé multi-chromosomique
    def extract_chr_ids(sequence: str) -> set[str]:
        ids = set()
        for part in sequence.split(';'):
            clean = re.sub(r"\D", "", part)
            if clean:
                ids.add(clean)
        return ids

    multi_der: dict[str, list[dict[str, object]]] = {}
    for an in norm_counts:
        if an.startswith(('der', 'dic')):
            m = re.match(r"^(?:der|dic)\(([0-9?;]+)\)", an)
            if m:
                # Inclut der(...), t(...), et les écritures déroulées
                # der(...)(8pter->...::17q11->17qter).
                chrs = {c for c in get_chromosomes(an) if c != '?'}
                if len(chrs) > 1:
                    centromeric = bool(re.search(r'[pq]10', an))
                    index = first_index.get(an, float('inf'))
                    for c in chrs:
                        multi_der.setdefault(c, []).append({
                            "anom": an,
                            "index": index,
                            "centromeric": centromeric,
                            "clones": norm_to_clones.get(an, set()),
                        })

    for an in norm_counts:
        if an.startswith(('+', '-')):
            num = re.sub(r"\D", "", an)
            if not num or num not in multi_der:
                continue
            entries = sorted(
                [
                    entry
                    for entry in multi_der[num]
                    if not clone_map
                    or norm_to_clones.get(an, set()).intersection(
                        entry.get("clones", set())  # type: ignore[arg-type]
                    )
                ],
                key=lambda e: e["index"],  # type: ignore[index]
            )
            if not entries:
                continue
            if an.startswith('-'):
                ref_entry = entries[0]
                reason = "Perte implicite"
            else:
                ref_entry = next(
                    (e for e in entries if e.get("centromeric")),  # type: ignore[arg-type]
                    None,
                )
                if not ref_entry:
                    continue
                reason = "Gain implicite (points de cassure p10/q10)"

            ref_norm = ref_entry["anom"]  # type: ignore[index]
            ref = norm_to_orig.get(ref_norm, ref_norm)
            implicit[an] = {"reason": reason, "ref": ref}

    # 3) Répétitions de la même anomalie (même chromosome)
    base_pattern = re.compile(
        r'^(?:der|dic|t|i|ider|idic|r)\([0-9;]+\)'
    )
    base_map = {}
    structural_refs: dict[tuple[str, str], list[tuple[str, str | None]]] = {}
    der_t_refs: dict[tuple[str, str, tuple[str, ...]], list[tuple[str, str | None]]] = {}

    def structural_repeat_key(anom: str) -> tuple[tuple[str, str], str | None] | None:
        m = re.match(r"^(del|dup|add|ins|inv)\(([^)]*)\)((?:\([^)]*\))*)$", anom, re.IGNORECASE)
        if not m:
            return None
        kind, chroms, breakpoints = m.groups()
        bp = breakpoints or None
        return (kind.lower(), chroms), bp

    def der_t_repeat_key(anom: str) -> tuple[tuple[str, str, tuple[str, ...]], str | None] | None:
        base = strip_multiplicity(anom)
        m = re.match(
            r"^(der|dic)\(([^)]*)\).*?t\(([^)]*)\)(?:\(([^)]*)\))?",
            base,
            re.IGNORECASE,
        )
        if not m:
            return None
        kind, anchor, t_chroms_raw, bp = m.groups()
        t_chroms = tuple(
            sorted((chrom.lstrip("?") or "?") for chrom in t_chroms_raw.split(";") if chrom)
        )
        return (kind.lower(), anchor.lstrip("?") or "?", t_chroms), bp

    def compatible_breakpoints(a: str | None, b: str | None) -> bool:
        return a is None or b is None or a == b

    for a in unique_anoms:
        norm = normalize_anomaly(a)
        # les gains/pertes répétés dans un même clone (ex: +8,+8)
        # correspondent à une trisomie ou tetrasomie et ne doivent pas
        # être considérés comme des duplications implicites
        if norm.startswith(('+', '-')):
            continue
        if norm.startswith(('der', 'dic')):
            parsed = der_t_repeat_key(norm)
            if parsed:
                key, bp = parsed
                refs = der_t_refs.setdefault(key, [])
                ref_norm = next(
                    (
                        ref_norm
                        for ref_norm, ref_bp in refs
                        if compatible_breakpoints(bp, ref_bp)
                    ),
                    None,
                )
                if ref_norm:
                    ref = norm_to_orig.get(ref_norm, ref_norm)
                    implicit.setdefault(
                        norm,
                        {"reason": "Duplication avec l'anomalie de référence", "ref": ref},
                    )
                else:
                    refs.append((norm, bp))
            # Les autres dérivés sont traités séparément via la comparaison des
            # ensembles chromosomiques. Ne pas les inclure ici pour éviter
            # de marquer implicites des dérivés distincts.
            continue
        if norm.startswith(("del", "dup", "add", "ins", "inv")):
            parsed = structural_repeat_key(norm)
            if not parsed:
                continue
            key, bp = parsed
            refs = structural_refs.setdefault(key, [])
            ref_norm = next(
                (
                    ref_norm
                    for ref_norm, ref_bp in refs
                    if compatible_breakpoints(bp, ref_bp)
                ),
                None,
            )
            if ref_norm:
                ref = norm_to_orig.get(ref_norm, ref_norm)
                implicit.setdefault(
                    norm,
                    {"reason": "Duplication avec l'anomalie de référence", "ref": ref},
                )
            else:
                refs.append((norm, bp))
            continue
        else:
            m = base_pattern.match(norm)
            base = m.group(0) if m else norm
        base_map.setdefault(base, []).append(norm)
    for norms in base_map.values():
        if len(norms) > 1:
            ref_norm = norms[0]
            ref = norm_to_orig.get(ref_norm, ref_norm)
            for n in norms[1:]:
                if n not in implicit:
                    implicit[n] = {"reason": "Duplication avec l'anomalie de référence", "ref": ref}

    # 4) Dérivés distincts impliquant le même ensemble de chromosomes
    derived_groups = {}
    for an in norm_counts:
        if an.startswith(('der', 'dic')):
            # Ne pas marquer implicites les dérivés porteurs d'une translocation
            # (cas des t déséquilibrées potentiellement légitimes).
            if "t(" in an:
                continue
            chroms = tuple(sorted(get_chromosomes(an)))
            if chroms:
                derived_groups.setdefault(chroms, []).append(an)

    for ders in derived_groups.values():
        if len(ders) <= 1:
            continue
        ders_sorted = sorted(ders, key=lambda d: first_index.get(d, float('inf')))
        ref_norm = ders_sorted[0]
        ref = norm_to_orig.get(ref_norm, ref_norm)
        for n in ders_sorted[1:]:
            implicit.setdefault(n, {"reason": "Dérivé implicite (mêmes chromosomes)", "ref": ref})

    return implicit


# =========================
# Scoring
# =========================
@dataclass(frozen=True)
class RuleSpec:
    rule_id: str
    system: str
    default_score: int
    title: str
    explanation: str


RULE_CATALOG_SHEET_ENV = "RULE_CATALOG_SHEET_URL"
RULE_EDITABLE_COLUMNS = ("Libellé", "Explication")
_RULE_TEXT_OVERRIDES_CACHE: dict[str, dict[str, dict[str, str]]] = {}


RULE_CATALOG: tuple[RuleSpec, ...] = (
    RuleSpec("ISCN.CONSTITUTIONAL_GAIN", "ISCN 2024", 0, "Gain constitutionnel", "Gain chromosomique annoté constitutionnel."),
    RuleSpec("ISCN.CONSTITUTIONAL_SUSPECT", "ISCN 2024", 0, "Suspicion constitutionnelle", "Anomalie annotée comme possiblement constitutionnelle."),
    RuleSpec("ISCN.CONSTITUTIONAL_CERTAIN", "ISCN 2024", 0, "Anomalie constitutionnelle certaine", "Anomalie annotée constitutionnelle."),
    RuleSpec("ISCN.REPEAT_NOTATION", "ISCN 2024", 0, "Notation de répétition", "Notation idem/sl/sdl déjà portée par un clone précédent."),
    RuleSpec("ISCN.ABSENT_IN_REPEAT_ZERO", "ISCN 2024", 0, "Absence dans un clone répété", "Anomalie indiquée absente dans un clone secondaire et déjà comptée auparavant."),
    RuleSpec("ISCN.ABSENT_IN_REPEAT_IMPLICIT_LOSS", "ISCN 2024", 1, "Perte implicite dans un clone répété", "Absence d'une anomalie structurale interprétée comme perte chromosomique implicite."),
    RuleSpec("ISCN.IMPLICIT", "ISCN 2024", 0, "Anomalie implicite", "Anomalie déjà expliquée par une anomalie de référence."),
    RuleSpec("ISCN.BALANCED_T", "ISCN 2024", 1, "Translocation équilibrée", "Translocation avec dérivés réciproques compatibles."),
    RuleSpec("ISCN.UNBALANCED_T", "ISCN 2024", 2, "Translocation déséquilibrée", "Translocation non équilibrée détectée explicitement."),
    RuleSpec("ISCN.T_ALREADY_IN_DER", "ISCN 2024", 0, "Translocation déjà comptée", "Translocation explicite déjà incluse dans un chromosome dérivé."),
    RuleSpec("ISCN.T_VIA_DER_BALANCED", "ISCN 2024", 1, "Translocation équilibrée via dérivé", "Translocation équilibrée comptée à partir d'un chromosome dérivé."),
    RuleSpec("ISCN.T_VIA_DER_UNBALANCED", "ISCN 2024", 2, "Translocation déséquilibrée via dérivé", "Translocation déséquilibrée comptée à partir d'un chromosome dérivé."),
    RuleSpec("ISCN.DER_T_COUNTED", "ISCN 2024", 0, "Dérivé associé déjà compté", "Dérivé associé à une translocation déjà comptabilisée."),
    RuleSpec("ISCN.DER_BALANCED_T", "ISCN 2024", 0, "Dérivé de translocation équilibrée", "Dérivé issu d'une translocation équilibrée déjà représentée."),
    RuleSpec("ISCN.INTERTWINED_DER_BALANCED_SIMPLE", "ISCN 2024", 1, "Dérivés enchevêtrés équilibrés", "Groupe de dérivés reliés avec cassures concordantes."),
    RuleSpec("ISCN.INTERTWINED_DER_BALANCED_INSERTION", "ISCN 2024", 2, "Dérivés enchevêtrés équilibrés avec insertion", "Groupe de dérivés reliés avec cassures concordantes et insertion additionnelle."),
    RuleSpec("ISCN.INTERTWINED_DER_UNBALANCED", "ISCN 2024", 2, "Dérivés enchevêtrés déséquilibrés", "Groupe de dérivés reliés avec cassures non concordantes."),
    RuleSpec("ISCN.INTERTWINED_DER_PART", "ISCN 2024", 0, "Partie d'un remaniement", "Dérivé déjà compté dans un remaniement enchevêtré."),
    RuleSpec("ISCN.STRUCTURAL_GAIN_DUPLICATE", "ISCN 2024", 1, "Gain structural dupliqué", "Gain d'une anomalie structurale déjà décrite."),
    RuleSpec("ISCN.SEMANTIC_PLUS_ISO", "ISCN 2024", 2, "Gain d'isochromosome", "Notation +i(...) interprétée comme gain chromosomique plus isochromosome."),
    RuleSpec("ISCN.SEMANTIC_PLUS_DEL", "ISCN 2024", 2, "Gain de chromosome délété", "Notation +del(...) interprétée comme gain chromosomique plus délétion."),
    RuleSpec("ISCN.MAR", "ISCN 2024", 1, "Marqueur", "Chromosome marqueur."),
    RuleSpec("ISCN.DICENTRIC", "ISCN 2024", 2, "Chromosome dicentrique", "Chromosome dicentrique compté comme remaniement complexe."),
    RuleSpec("ISCN.DER_NO_BREAKPOINT", "ISCN 2024", 1, "Dérivé sans cassure", "Chromosome dérivé sans point de cassure détaillé."),
    RuleSpec("ISCN.DER_MULTI", "ISCN 2024", 2, "Dérivé multichromosomique", "Chromosome dérivé impliquant au moins deux chromosomes identifiés."),
    RuleSpec("ISCN.DER_MULTI_UNCERTAIN", "ISCN 2024", 2, "Dérivé multichromosomique incertain", "Dérivé impliquant plusieurs chromosomes avec positions incertaines."),
    RuleSpec("ISCN.DER_UNCERTAIN_SECOND", "ISCN 2024", 1, "Dérivé avec second chromosome incertain", "Dérivé sans certitude sur l'implication d'un second chromosome."),
    RuleSpec("ISCN.DER_SAME_CHR", "ISCN 2024", 1, "Dérivé intrachromosomique", "Chromosome dérivé issu du même chromosome."),
    RuleSpec("ISCN.SINGLE_CHR_GAIN_REPEAT", "ISCN 2024", 2, "Gain chromosomique répété", "Plusieurs gains du même chromosome dans la formule."),
    RuleSpec("ISCN.SINGLE_CHR_TRIPLICATION", "ISCN 2024", 2, "Triplication", "Anomalie de triplication d'un chromosome ou segment."),
    RuleSpec("ISCN.SINGLE_CHR_ISODERIVATIVE", "ISCN 2024", 2, "Isodérivé", "Chromosome isodérivé ou isodicentrique."),
    RuleSpec("ISCN.COMPLEX_MULTI_CHR", "ISCN 2024", 2, "Déséquilibre multichromosomique", "Anomalie déséquilibrée impliquant plusieurs chromosomes."),
    RuleSpec("ISCN.UNBALANCED_TRANSLOCATION", "ISCN 2024", 2, "Translocation déséquilibrée", "Translocation non pure ou portée par un dérivé."),
    RuleSpec("ISCN.GAIN_SIMPLE", "ISCN 2024", 1, "Gain simple", "Gain chromosomique simple."),
    RuleSpec("ISCN.LOSS_SIMPLE", "ISCN 2024", 1, "Perte simple", "Perte chromosomique simple."),
    RuleSpec("ISCN.OTHER_STANDARD", "ISCN 2024", 1, "Autre anomalie standard", "Anomalie non constitutionnelle non couverte par une règle plus spécifique."),
    RuleSpec("JON.CONSTITUTIONAL_GAIN", "Jondreville 2020", 0, "Gain constitutionnel", "Gain chromosomique annoté constitutionnel."),
    RuleSpec("JON.CONSTITUTIONAL_SUSPECT", "Jondreville 2020", 0, "Suspicion constitutionnelle", "Anomalie annotée comme possiblement constitutionnelle."),
    RuleSpec("JON.CONSTITUTIONAL_CERTAIN", "Jondreville 2020", 0, "Anomalie constitutionnelle certaine", "Anomalie annotée constitutionnelle."),
    RuleSpec("JON.REPEAT_NOTATION", "Jondreville 2020", 0, "Notation de répétition", "Notation idem/sl/sdl déjà portée par un clone précédent."),
    RuleSpec("JON.ABSENT_IN_REPEAT_ZERO", "Jondreville 2020", 0, "Absence dans un clone répété", "Anomalie indiquée absente dans un clone secondaire et déjà comptée auparavant."),
    RuleSpec("JON.ABSENT_IN_REPEAT_IMPLICIT_LOSS", "Jondreville 2020", 1, "Perte implicite dans un clone répété", "Absence d'une anomalie structurale interprétée comme perte chromosomique implicite."),
    RuleSpec("JON.IMPLICIT", "Jondreville 2020", 0, "Anomalie implicite", "Duplication avec une anomalie de référence."),
    RuleSpec("JON.MAR", "Jondreville 2020", 1, "Marqueur", "Chromosome marqueur."),
    RuleSpec("JON.TRIPLOIDY", "Jondreville 2020", 0, "Triploïdie", "Triploïdie ignorée dans le calcul Jondreville."),
    RuleSpec("JON.DEFAULT", "Jondreville 2020", 1, "Anomalie non constitutionnelle", "Chaque anomalie non constitutionnelle vaut un point."),
)

RULE_PRIORITY: dict[str, tuple[str, ...]] = {
    "ISCN 2024": (
        "ISCN.ABSENT_IN_REPEAT_ZERO",
        "ISCN.ABSENT_IN_REPEAT_IMPLICIT_LOSS",
        "ISCN.CONSTITUTIONAL_GAIN",
        "ISCN.CONSTITUTIONAL_SUSPECT",
        "ISCN.CONSTITUTIONAL_CERTAIN",
        "ISCN.REPEAT_NOTATION",
        "ISCN.INTERTWINED_DER_BALANCED_SIMPLE",
        "ISCN.INTERTWINED_DER_BALANCED_INSERTION",
        "ISCN.INTERTWINED_DER_UNBALANCED",
        "ISCN.INTERTWINED_DER_PART",
        "ISCN.BALANCED_T",
        "ISCN.UNBALANCED_T",
        "ISCN.T_ALREADY_IN_DER",
        "ISCN.IMPLICIT",
        "ISCN.STRUCTURAL_GAIN_DUPLICATE",
        "ISCN.SEMANTIC_PLUS_ISO",
        "ISCN.SEMANTIC_PLUS_DEL",
        "ISCN.MAR",
        "ISCN.DICENTRIC",
        "ISCN.T_VIA_DER_BALANCED",
        "ISCN.T_VIA_DER_UNBALANCED",
        "ISCN.DER_T_COUNTED",
        "ISCN.DER_BALANCED_T",
        "ISCN.DER_NO_BREAKPOINT",
        "ISCN.DER_MULTI",
        "ISCN.DER_MULTI_UNCERTAIN",
        "ISCN.DER_UNCERTAIN_SECOND",
        "ISCN.DER_SAME_CHR",
        "ISCN.SINGLE_CHR_GAIN_REPEAT",
        "ISCN.SINGLE_CHR_TRIPLICATION",
        "ISCN.SINGLE_CHR_ISODERIVATIVE",
        "ISCN.COMPLEX_MULTI_CHR",
        "ISCN.UNBALANCED_TRANSLOCATION",
        "ISCN.GAIN_SIMPLE",
        "ISCN.LOSS_SIMPLE",
        "ISCN.OTHER_STANDARD",
    ),
    "Jondreville 2020": (
        "JON.ABSENT_IN_REPEAT_ZERO",
        "JON.ABSENT_IN_REPEAT_IMPLICIT_LOSS",
        "JON.CONSTITUTIONAL_GAIN",
        "JON.CONSTITUTIONAL_SUSPECT",
        "JON.CONSTITUTIONAL_CERTAIN",
        "JON.REPEAT_NOTATION",
        "JON.IMPLICIT",
        "JON.MAR",
        "JON.TRIPLOIDY",
        "JON.DEFAULT",
    ),
}

RULE_BY_ID = {rule.rule_id: rule for rule in RULE_CATALOG}

RULE_TECHNICAL_CHECKS: dict[str, str] = {
    "ISCN.CONSTITUTIONAL_GAIN": "Le code normalise l'anomalie puis teste la regex ^\\+\\d+c$. En langage courant: le texte doit commencer par '+', contenir ensuite un ou plusieurs chiffres, puis finir par 'c'. Exemple: '+8c'. Le 'c' indique constitutionnel. Cette règle concerne donc un gain chromosomique annoté constitutionnel. Score 0.",
    "ISCN.CONSTITUTIONAL_SUSPECT": "Le code normalise l'anomalie puis vérifie si le texte finit exactement par '?c'. Le '?' indique une incertitude dans la notation et le 'c' indique constitutionnel; ensemble, le code l'interprète comme suspicion d'anomalie constitutionnelle. Cette règle est séparée pour permettre une explication clinique propre. Score 0.",
    "ISCN.CONSTITUTIONAL_CERTAIN": "Le code normalise l'anomalie puis vérifie si le texte finit par 'c', après avoir déjà exclu les cas '+<nombre>c' et '?c'. Cela correspond aux autres anomalies annotées constitutionnelles certaines. Score 0.",
    "ISCN.REPEAT_NOTATION": "Le code met le texte en minuscules puis teste les notations qui ne décrivent pas une nouvelle anomalie mais renvoient à un clone précédent. 'idem' est accepté tel quel. La regex ^(?:sl|sdl)\\d*$ veut dire: début du texte (^), puis soit 'sl' soit 'sdl' ((?:sl|sdl)), puis éventuellement des chiffres (\\d*), puis fin du texte ($). Exemples retenus: sl, sl2, sdl, sdl3. Score 0 car l'anomalie est déjà portée par une notation précédente.",
    "ISCN.ABSENT_IN_REPEAT_ZERO": "Avant le scoring, l'analyse des clones repère une anomalie héritée d'un clone précédent mais explicitement absente dans le clone courant. Elle est stockée dans zeroed_reasons avec score_override=0. Le score est donc forcé à 0 et les règles suivantes ne sont pas testées pour cette anomalie.",
    "ISCN.ABSENT_IN_REPEAT_IMPLICIT_LOSS": "Même mécanisme que ABSENT_IN_REPEAT_ZERO, mais l'anomalie absente est une anomalie structurale dont l'absence est interprétée comme perte chromosomique implicite: der(...), dic(...), add(...) ou certains r(...). Elle est stockée avec score_override=1. Score 1.",
    "ISCN.IMPLICIT": "Le code cherche si l'anomalie est déjà expliquée par une autre anomalie de référence. Il normalise les écritures, compare les répétitions structurales del/dup/add/ins/inv avec points de cassure compatibles, repère les gains/pertes simples expliqués par un dérivé multi-chromosomique, et regroupe certains dérivés qui impliquent le même ensemble de chromosomes. Si une référence existe, cette anomalie vaut 0 pour éviter un double comptage.",
    "ISCN.BALANCED_T": "Le code traite d'abord une anomalie écrite t(...). Il extrait les chromosomes de la translocation, par exemple t(9;22) donne la clé {9,22}. Ensuite il cherche dans le même clone des der(...)/ider(...) qui portent la même translocation. La règle équilibrée est retenue seulement s'il existe au moins deux chromosomes d'ancrage différents, au moins un der(...), et des points de cassure absents ou identiques, sans '?' dans les cassures. Score 1.",
    "ISCN.UNBALANCED_T": "Le code part aussi d'une anomalie t(...). Il extrait la même clé chromosomique que pour BALANCED_T. Si une relation avec un dérivé du même clone existe, mais que les critères d'équilibre ne sont pas remplis (ancrages insuffisants, cassures discordantes ou incertaines), la translocation est considérée déséquilibrée. Score 2.",
    "ISCN.T_ALREADY_IN_DER": "Quand une translocation t(...) est écrite explicitement, le code vérifie si un der(...) du même clone contient déjà une t(...) avec exactement les mêmes chromosomes. Exemple: t(9;22) et der(22)t(9;22) dans le même clone. Si oui, la translocation explicite vaut 0 car elle est déjà représentée par le dérivé.",
    "ISCN.T_VIA_DER_BALANCED": "Si l'anomalie est un der(...) contenant une t(...), le code extrait les chromosomes de cette t(...), vérifie qu'au moins deux chromosomes connus sont présents, sans '?', et que la translocation n'a pas déjà été comptée. Si la même clé de translocation est reconnue comme équilibrée dans le clone, le dérivé sert à compter une translocation équilibrée. Score 1.",
    "ISCN.T_VIA_DER_UNBALANCED": "Même détection que T_VIA_DER_BALANCED, mais la clé de translocation portée par le dérivé n'est pas reconnue comme équilibrée dans le clone. Le dérivé sert donc à compter une translocation déséquilibrée. Score 2.",
    "ISCN.DER_T_COUNTED": "Pour un der(...) contenant une t(...), le code regarde si la même translocation a déjà été comptée dans le clone. Il compare une clé composée des chromosomes de t(...), par exemple {9,22}. Si cette clé est déjà présente dans les translocations explicites ou déjà comptées, le dérivé vaut 0 pour éviter de compter deux fois le même événement.",
    "ISCN.DER_BALANCED_T": "Le code identifie un der(...) qui fait partie d'une translocation équilibrée déjà reconstruite par paire de dérivés. Il vérifie aussi que le dérivé ne contient pas de composant additionnel comme add(...), del(...), dup(...), ins(...), inv(...) ou r(...). Si c'est juste la représentation du dérivé équilibré déjà pris en compte, score 0.",
    "ISCN.INTERTWINED_DER_BALANCED_SIMPLE": "Dans un même clone, le code regroupe au moins trois der(...)/dic(...) qui partagent des chromosomes. Il exige au moins trois chromosomes impliqués, trois ancrages cohérents, et des signatures 'chromosome + point de cassure' retrouvées au moins deux fois et couvrant le groupe. Sans insertion détectée, le groupe équilibré vaut 1.",
    "ISCN.INTERTWINED_DER_BALANCED_INSERTION": "Même groupe équilibré que INTERTWINED_DER_BALANCED_SIMPLE, mais le code détecte aussi une insertion: présence de ins(...) ou implication d'au moins trois chromosomes dans un composant. Le score est augmenté à 2 pour porter cette complexité additionnelle.",
    "ISCN.INTERTWINED_DER_UNBALANCED": "Même recherche de groupe que INTERTWINED_DER_BALANCED: au moins trois der(...)/dic(...) reliés dans le même clone. Ici, les cassures ne sont pas suffisamment concordantes pour conclure à un équilibre exact. Le code accepte au moins deux signatures de cassure supportées, y compris avec incertitude '?'. Le groupe est alors considéré déséquilibré. Score 2.",
    "ISCN.INTERTWINED_DER_PART": "Une fois qu'un groupe de dérivés enchevêtrés est détecté, le code choisit la première anomalie du groupe comme porteuse du score. Les autres dérivés du même groupe reçoivent cette règle avec score 0, car ils font partie du même remaniement déjà compté.",
    "ISCN.STRUCTURAL_GAIN_DUPLICATE": "Le code ne regarde cette règle que si l'anomalie commence par '+'. Il enlève le '+' et la multiplicité éventuelle xN/×N, puis obtient une base structurale comme der(...), dic(...), del(...), dup(...), inv(...), ins(...), add(...), i(...), ider(...), idic(...) ou r(...). Il cherche si la même base structurale est déjà apparue avant dans la formule. Si oui, le '+' correspond à un gain d'une anomalie structurale déjà décrite. Score 1.",
    "ISCN.SEMANTIC_PLUS_ISO": "Le code traite une anomalie commençant par '+'. Après retrait du '+', si la base commence par i(...), il interprète '+i(...)' comme deux informations: gain du chromosome concerné plus isochromosome. Score 2.",
    "ISCN.SEMANTIC_PLUS_DEL": "Le code traite une anomalie commençant par '+'. Après retrait du '+', si la base commence par del(...), il interprète '+del(...)' comme deux informations: gain du chromosome concerné plus délétion. Score 2.",
    "ISCN.MAR": "Le code enlève d'abord le signe '+' ou '-' et une multiplicité finale xN/×N. Il met en minuscules puis teste si la base commence par 'mar' ou finit par 'mar'. Il reconnaît aussi les quantités implicites: '+2mar' signifie deux marqueurs, '+1~3mar' prend le maximum 3, '+mar1x2' signifie le marqueur 1 répété deux fois. Chaque marqueur retenu est scoré selon la logique marqueur.",
    "ISCN.DICENTRIC": "Après normalisation, retrait du signe et retrait d'une multiplicité xN/×N, le code teste simplement si le texte commence par 'dic'. Cela correspond à une notation de chromosome dicentrique, par exemple dic(9;20). Score 2.",
    "ISCN.DER_NO_BREAKPOINT": "Le code applique la regex ^der\\([^)]+\\)$ sur la base sans signe ni multiplicité. Elle veut dire: début du texte (^), puis 'der(', puis au moins un caractère qui n'est pas ')' ([^)]+), puis ')', puis fin du texte ($). Donc der(7) passe, mais der(7)(p10q10) ne passe pas car il contient une seconde parenthèse de cassure. Score 1.",
    "ISCN.DER_MULTI": "Pour un der(...), le code extrait les chromosomes mentionnés après des mots-clés comme der, dic, del, dup, ins, t, i, ider, idic ou r, et aussi les chromosomes écrits dans la seconde parenthèse des dérivés. Les '?' sont exclus du compte des chromosomes connus. Si au moins deux chromosomes connus sont trouvés et que la notation n'est pas incertaine, le dérivé est multi-chromosomique. Score 2.",
    "ISCN.DER_MULTI_UNCERTAIN": "Même extraction que DER_MULTI, mais le texte contient un '?'. Le code garde la règle multi-chromosomique si au moins deux chromosomes connus sont présents, tout en signalant que la notation est incertaine. Score 2, avec une explication d'imprécision.",
    "ISCN.DER_UNCERTAIN_SECOND": "Pour un der(...), le code extrait les chromosomes connus et cherche le caractère '?'. Si un seul chromosome certain est retrouvé mais qu'un '?' est présent, le code comprend qu'un second chromosome pourrait être impliqué sans certitude. Score 1.",
    "ISCN.DER_SAME_CHR": "Cette règle est le cas restant pour der(...). Elle arrive après exclusion des translocations déjà traitées, des der(...) sans cassure, des dérivés multi-chromosomiques et des dérivés avec second chromosome incertain. Le code considère alors que le dérivé concerne le même chromosome ou un remaniement intrachromosomique. Score 1.",
    "ISCN.SINGLE_CHR_GAIN_REPEAT": "Le code détecte une anomalie commençant par '+' répétée plus d'une fois après normalisation, par exemple deux occurrences de +8. Cela correspond à un gain répété du même chromosome. Score 2, plafonné pour éviter de multiplier artificiellement le score par le nombre d'occurrences.",
    "ISCN.SINGLE_CHR_TRIPLICATION": "Après retrait éventuel du signe et normalisation, le code teste si l'anomalie commence par 'trp'. Cela correspond à une triplication. Score 2.",
    "ISCN.SINGLE_CHR_ISODERIVATIVE": "Après normalisation, le code teste si l'anomalie commence par 'ider'. Cela correspond à un chromosome isodérivé ou isodicentrique. Score 2.",
    "ISCN.COMPLEX_MULTI_CHR": "Le code extrait d'abord les chromosomes connus. S'il y en a zéro ou un, la règle ne s'applique pas. S'il y en a au moins deux, il teste ensuite des formes déséquilibrées complexes: dic(...), r(...), insertion non équilibrée, ou translocation non équilibrée. Pour der(...), il exclut d'abord les formes der(...) sans cassure car elles ne prouvent pas clairement plusieurs chromosomes. Score 2.",
    "ISCN.UNBALANCED_TRANSLOCATION": "Le code considère deux situations. Première situation: le texte contient der(...) ou dic(...) et contient aussi t(...), donc une translocation portée par un dérivé/dicentrique. Deuxième situation: le texte contient t(...) mais n'est pas une translocation pure équilibrée de forme t(chromosome;chromosome)(points de cassure), sans der, sans '+' et sans '-'. Score 2.",
    "ISCN.GAIN_SIMPLE": "Cette règle arrive tard. Si aucune règle plus spécifique n'a retenu l'anomalie et que le texte normalisé commence par '+', le code la considère comme gain chromosomique simple. Score 1.",
    "ISCN.LOSS_SIMPLE": "Cette règle arrive tard. Si aucune règle plus spécifique n'a retenu l'anomalie et que le texte normalisé commence par '-', le code la considère comme perte chromosomique simple. Score 1.",
    "ISCN.OTHER_STANDARD": "Dernier recours ISCN. Si l'anomalie n'est pas constitutionnelle, pas répétée, pas implicite, pas marqueur, pas dérivé traité, pas déséquilibre identifié par les règles précédentes, le code applique le score standard. Score 1.",
    "JON.CONSTITUTIONAL_GAIN": "Même détection que ISCN.CONSTITUTIONAL_GAIN: regex ^\\+\\d+c$, c'est-à-dire '+', puis un ou plusieurs chiffres, puis 'c' en fin de texte. Exemple '+8c'. Jondreville attribue 0 point.",
    "JON.CONSTITUTIONAL_SUSPECT": "Même détection que ISCN.CONSTITUTIONAL_SUSPECT: le texte normalisé finit par '?c'. Le code distingue ce cas pour afficher une explication clinique spécifique à une suspicion constitutionnelle. Score 0.",
    "JON.CONSTITUTIONAL_CERTAIN": "Même détection que ISCN.CONSTITUTIONAL_CERTAIN: le texte normalisé finit par 'c', hors cas '+<nombre>c' et '?c' déjà traités. Score 0.",
    "JON.REPEAT_NOTATION": "Même logique que ISCN.REPEAT_NOTATION. Le code reconnaît 'idem' ou les formes sl/sdl avec numéro optionnel par la regex ^(?:sl|sdl)\\d*$. Ce sont des renvois à un clone précédent, pas de nouvelles anomalies à scorer. Score 0.",
    "JON.ABSENT_IN_REPEAT_ZERO": "Même source que ISCN.ABSENT_IN_REPEAT_ZERO: anomalie héritée d'un clone précédent puis indiquée absente dans le clone courant. zeroed_reasons impose score_override=0. Score 0.",
    "JON.ABSENT_IN_REPEAT_IMPLICIT_LOSS": "Même source que ISCN.ABSENT_IN_REPEAT_IMPLICIT_LOSS: absence d'une anomalie structurale interprétée comme perte chromosomique implicite dans un clone répété. zeroed_reasons impose score_override=1. Score 1.",
    "JON.IMPLICIT": "Pour Jondreville, le code ne met à 0 que les anomalies classées comme 'Duplication avec l'anomalie de référence'. Cela évite de recompter une anomalie déjà représentée dans un clone ou une anomalie structurale précédente. Les autres implicites ne sont pas tous traités de la même façon qu'en ISCN.",
    "JON.MAR": "Même détection technique que ISCN.MAR: retrait du signe et de la multiplicité, passage en minuscules, puis test 'commence par mar' ou 'finit par mar'. Les marqueurs sont reconnus même dans des formes quantifiées comme +2mar ou +mar1x2. Score 1 par règle Jondreville.",
    "JON.TRIPLOIDY": "Après normalisation complète, le code compare exactement le texte à 'triploidy' en minuscules. Il n'utilise pas de synonymes ni de regex ici: seule cette chaîne exacte déclenche la règle. Score 0 dans le calcul Jondreville.",
    "JON.DEFAULT": "Dernier recours Jondreville. Si aucune règle précédente ne s'applique, donc pas constitutionnel, pas répétition, pas implicite retenu, pas marqueur et pas triploïdie, le code compte l'anomalie comme anomalie non constitutionnelle standard. Score 1.",
}


def _google_sheet_csv_url(url_or_id: str) -> str:
    """Construit l'URL CSV publique d'un Google Sheet depuis une URL ou un ID."""

    source = str(url_or_id or "").strip()
    if not source:
        return ""
    if "docs.google.com" not in source:
        return f"https://docs.google.com/spreadsheets/d/{source}/export?format=csv&gid=0"

    parsed = urlparse(source)
    match = re.search(r"/d/([A-Za-z0-9_-]+)", parsed.path)
    if not match:
        return source
    sheet_id = match.group(1)
    query = parse_qs(parsed.query)
    gid = query.get("gid", ["0"])[0]
    if not gid and parsed.fragment:
        gid = parse_qs(parsed.fragment).get("gid", ["0"])[0]
    return f"https://docs.google.com/spreadsheets/d/{sheet_id}/export?format=csv&gid={gid or '0'}"


def _normalize_rule_catalog_columns(df: pd.DataFrame) -> pd.DataFrame:
    aliases = {
        "rule id": "Rule ID",
        "rule_id": "Rule ID",
        "id": "Rule ID",
        "libelle": "Libellé",
        "libellé": "Libellé",
        "label": "Libellé",
        "title": "Libellé",
        "explication": "Explication",
        "explanation": "Explication",
    }
    return df.rename(
        columns={
            col: aliases.get(str(col).strip().lower(), str(col).strip())
            for col in df.columns
        }
    )


def load_rule_catalog_public_texts(url_or_id: str | None = None) -> dict[str, dict[str, str]]:
    """Charge uniquement les textes publics modifiables du catalogue."""

    source = url_or_id or os.environ.get(RULE_CATALOG_SHEET_ENV, "")
    if not source:
        return {}
    if source in _RULE_TEXT_OVERRIDES_CACHE:
        return _RULE_TEXT_OVERRIDES_CACHE[source]

    df = pd.read_csv(_google_sheet_csv_url(source))
    df = _normalize_rule_catalog_columns(df)
    if "Rule ID" not in df.columns:
        raise ValueError("La feuille des règles doit contenir une colonne 'Rule ID'.")

    allowed_ids = set(RULE_BY_ID)
    overrides: dict[str, dict[str, str]] = {}
    unknown_ids: list[str] = []
    for _, row in df.iterrows():
        rule_id = str(row.get("Rule ID", "")).strip()
        if not rule_id:
            continue
        if rule_id not in allowed_ids:
            unknown_ids.append(rule_id)
            continue

        editable: dict[str, str] = {}
        for col in RULE_EDITABLE_COLUMNS:
            if col not in df.columns:
                continue
            value = row.get(col)
            if pd.isna(value):
                continue
            text = str(value).strip()
            if text:
                editable[col] = text
        if editable:
            overrides[rule_id] = editable

    if unknown_ids:
        raise ValueError(
            "La feuille des règles contient des Rule ID inconnus: "
            + ", ".join(sorted(set(unknown_ids)))
        )
    _RULE_TEXT_OVERRIDES_CACHE[source] = overrides
    return overrides


def _public_rule_text_overrides() -> dict[str, dict[str, str]]:
    try:
        return load_rule_catalog_public_texts()
    except Exception:
        return {}


def _rule_with_public_text(rule: RuleSpec, overrides: dict[str, dict[str, str]]) -> RuleSpec:
    override = overrides.get(rule.rule_id, {})
    return RuleSpec(
        rule.rule_id,
        rule.system,
        rule.default_score,
        override.get("Libellé", rule.title),
        override.get("Explication", rule.explanation),
    )


def validate_rule_catalog_integrity() -> None:
    """Vérifie que le catalogue canonique et l'ordre de priorité restent cohérents."""

    catalog_ids = [rule.rule_id for rule in RULE_CATALOG]
    duplicate_ids = sorted({rule_id for rule_id in catalog_ids if catalog_ids.count(rule_id) > 1})
    if duplicate_ids:
        raise ValueError("Rule ID dupliqués: " + ", ".join(duplicate_ids))

    catalog_id_set = set(catalog_ids)
    priority_ids = {rule_id for priority in RULE_PRIORITY.values() for rule_id in priority}
    missing_from_catalog = sorted(priority_ids - catalog_id_set)
    if missing_from_catalog:
        raise ValueError(
            "Rule ID présents dans RULE_PRIORITY mais absents du catalogue: "
            + ", ".join(missing_from_catalog)
        )

    missing_from_priority = sorted(catalog_id_set - priority_ids)
    if missing_from_priority:
        raise ValueError(
            "Rule ID présents dans le catalogue mais absents de RULE_PRIORITY: "
            + ", ".join(missing_from_priority)
        )


def get_rule_catalog_dataframe(public_text_url: str | None = None) -> pd.DataFrame:
    """Retourne le référentiel des règles affichable/exportable."""

    overrides = (
        load_rule_catalog_public_texts(public_text_url)
        if public_text_url
        else _public_rule_text_overrides()
    )
    return pd.DataFrame(
        [
            {
                "Rule ID": displayed.rule_id,
                "Référentiel": displayed.system,
                "Score par défaut": displayed.default_score,
                "Critère technique": RULE_TECHNICAL_CHECKS.get(displayed.rule_id, ""),
                "Libellé": displayed.title,
                "Explication": displayed.explanation,
            }
            for displayed in (
                _rule_with_public_text(rule, overrides) for rule in RULE_CATALOG
            )
        ]
    )


def get_rule_path(rule_id: str, public_text_url: str | None = None) -> list[dict[str, object]]:
    """Retourne le parcours de priorité menant à une règle retenue."""

    if not rule_id:
        return []

    overrides = (
        load_rule_catalog_public_texts(public_text_url)
        if public_text_url
        else _public_rule_text_overrides()
    )
    rule = RULE_BY_ID.get(rule_id)
    if rule:
        rule = _rule_with_public_text(rule, overrides)
    system = rule.system if rule else ("Jondreville 2020" if rule_id.startswith("JON.") else "ISCN 2024")
    priority = list(RULE_PRIORITY.get(system, ()))
    if rule_id not in priority:
        priority.append(rule_id)

    path = []
    for index, candidate_id in enumerate(priority, start=1):
        candidate = RULE_BY_ID.get(candidate_id)
        if not candidate:
            title = candidate_id
            score = ""
            explanation = ""
        else:
            candidate = _rule_with_public_text(candidate, overrides)
            title = candidate.title
            score = candidate.default_score
            explanation = candidate.explanation
        selected = candidate_id == rule_id
        path.append(
            {
                "order": index,
                "rule_id": candidate_id,
                "title": title,
                "default_score": score,
                "technical_check": RULE_TECHNICAL_CHECKS.get(candidate_id, ""),
                "explanation": explanation,
                "selected": selected,
            }
        )
        if selected:
            break

    return path


@dataclass(frozen=True)
class RuleDecision:
    rule_id: str
    score: int
    explanation: str


def apply_rule(condition: bool, rule_id: str, score: int, explanation: str) -> RuleDecision | None:
    if condition:
        return RuleDecision(rule_id=rule_id, score=score, explanation=explanation)
    return None


def absent_in_repeat_rule_decision(prefix: str, score: int, reason: str) -> RuleDecision:
    rule_suffix = (
        "ABSENT_IN_REPEAT_IMPLICIT_LOSS"
        if score == 1 or "Perte chromosomique implicite" in reason
        else "ABSENT_IN_REPEAT_ZERO"
    )
    return RuleDecision(
        rule_id=f"{prefix}.{rule_suffix}",
        score=score,
        explanation=reason,
    )


def single_chr_deseq_rule_decision(norm: str, count: int) -> RuleDecision | None:
    if norm.startswith("+") and count > 1:
        return RuleDecision(
            rule_id="ISCN.SINGLE_CHR_GAIN_REPEAT",
            score=2,
            explanation="Gain chromosomique répété",
        )
    if norm.startswith("trp"):
        return RuleDecision(
            rule_id="ISCN.SINGLE_CHR_TRIPLICATION",
            score=2,
            explanation="Triplication",
        )
    if norm.startswith("ider"):
        return RuleDecision(
            rule_id="ISCN.SINGLE_CHR_ISODERIVATIVE",
            score=2,
            explanation="Isodérivé",
        )
    return None


def calcul_score_jondroville(anomalies, clone_map, entries=None, zeroed_reasons: dict[str, tuple[int, str]] | None = None):
    """Calcule le score global selon Jondreville 2020."""

    def jon_dedup_key(anom: str) -> str:
        norm = normalize_anomaly(anom)
        if norm.startswith("t("):
            m = re.match(r"^t\(([^)]+)\)", norm, re.IGNORECASE)
            if m:
                return f"t({m.group(1)})"
        return norm

    filtered = []
    jon_dup_expl: dict[str, str] = {}
    if entries is not None:
        seen_by_key: dict[str, set[str]] = {}
        for entry in entries:
            anom = entry["anomaly"]
            clone = entry["clone"]
            key = jon_dedup_key(anom)
            clones = seen_by_key.setdefault(key, set())
            # Déduplication uniquement entre clones différents
            if clones and clone not in clones:
                jon_dup_expl[anom] = "Duplication avec l'anomalie de référence"
                continue
            clones.add(clone)
            filtered.extend([anom] * int(entry.get("count", 1)))
    else:
        filtered = list(anomalies)

    counts = Counter(filtered)
    implicit_info = detect_implicit_anomalies(filtered)
    total = 0
    scores = {}
    explanations = {}

    rule_ids = {}
    rule_expls = {}

    for anom, cnt in counts.items():
        norm = normalize_anomaly(anom)
        eff_cnt = effective_occurrences(anom, cnt, clone_map)
        # Ignorer les anomalies constitutionnelles (+Nc)
        zeroed_entry = None
        if zeroed_reasons:
            zeroed_entry = zeroed_reasons.get(anom) or zeroed_reasons.get(normalize_anomaly(anom))
        if zeroed_entry:
            override_score, reason = zeroed_entry
            decision = absent_in_repeat_rule_decision("JON", override_score, reason)
        else:
            decision = constitutional_rule_decision("JON", norm)
        if decision is None:
            decision = apply_rule(
                is_repeat_notation(norm),
                "JON.REPEAT_NOTATION",
                0,
                "Anomalies déjà connues dans un autre clone",
            )
        if (
            decision is None
            and norm in implicit_info
            and implicit_info[norm]["reason"] == "Duplication avec l'anomalie de référence"
        ):
            info = implicit_info[norm]
            decision = RuleDecision(
                rule_id="JON.IMPLICIT",
                score=0,
                explanation=f"{info['reason']} ({info['ref']})",
            )
        if decision is None:
            decision = apply_rule(
                is_marker_anomaly(norm),
                "JON.MAR",
                1,
                "Marqueur",
            )
        if decision is None:
            decision = apply_rule(
                norm.lower() == 'triploidy',
                "JON.TRIPLOIDY",
                0,
                "Triploïdie ignorée dans le calcul",
            )
        if decision is None:
            decision = RuleDecision(
                rule_id="JON.DEFAULT",
                score=1,  # Chaque anomalie non-constitutionnelle vaut 1 point
                explanation="-",
            )

        score_per_occurrence = decision.score
        explanation = decision.explanation
        if anom in jon_dup_expl:
            explanation = jon_dup_expl[anom]
        explanation = append_uncertainty_note(anom, explanation)
        score = score_per_occurrence * eff_cnt
        if decision.rule_id.startswith("JON.ABSENT_IN_REPEAT"):
            if decision.score == 1:
                score = 1
            elif "Perte chromosomique implicite" in explanation:
                score = 1
        scores[anom] = score
        explanations[anom] = explanation
        total += score
        rule_ids[anom] = decision.rule_id
        rule_expls[anom] = decision.explanation

    return scores, explanations, total, rule_ids, rule_expls


def calcul_score_iscn(
    anomalies,
    clone_map,
    zeroed_reasons: dict[str, tuple[int, str]] | None = None,
    scoring_clone_map=None,
):
    """Calcule le détail des scores selon la grille ISCN 2024."""

    counts = Counter(anomalies)
    norm_counts = Counter(normalize_anomaly(a) for a in anomalies)
    implicit_info = detect_implicit_anomalies(
        anomalies, scoring_clone_map or clone_map
    )
    first_index = {}
    for idx, anom in enumerate(anomalies):
        first_index.setdefault(normalize_anomaly(anom), idx)

    rows = []
    total = 0

    rule_ids = {}
    rule_expls = {}

    # Détection des translocations équilibrées sous forme de dérivés "pairés"
    def extract_t_entries(anom_str: str):
        entries = []
        for m in re.finditer(r"t\(([^)]*)\)(?:\(([^)]*)\))?", anom_str):
            chroms_raw = m.group(1)
            bp_raw = m.group(2)
            chroms = [c.lstrip("?") for c in chroms_raw.split(";") if c]
            if not chroms:
                continue
            key = tuple(sorted(chroms))
            t_str = f"t({chroms_raw})"
            if bp_raw:
                t_str = f"{t_str}({bp_raw})"
            entries.append({"key": key, "bp": bp_raw, "t_str": t_str})
        return entries

    def has_extra_der_components(anom_str: str) -> bool:
        return bool(re.search(r"(add|del|dup|ins|inv|r\()", anom_str))

    def t_entries_in_str(anom_str: str) -> list[tuple[tuple[str, ...], str | None]]:
        entries: list[tuple[tuple[str, ...], str | None]] = []
        for m in re.finditer(r"t\(([^)]*)\)(?:\(([^)]*)\))?", anom_str, re.IGNORECASE):
            chroms = [c.lstrip("?") for c in m.group(1).split(";") if c]
            if not chroms:
                continue
            key = tuple(sorted(chroms))
            bp = m.group(2)
            entries.append((key, bp))
        return entries

    der_t_map: dict[tuple[str, ...], dict[str, object]] = {}
    der_t_by_anom: dict[str, str] = {}
    der_t_any: dict[str, set[tuple[str, ...]]] = {}
    ider_t_any: dict[str, set[tuple[str, ...]]] = {}
    # clone -> tkey -> list of (anom, anchor, is_der)
    t_pairs_by_clone: dict[str, dict[tuple[str, ...], list[tuple[str, str, bool]]]] = {}
    # clone -> tkey -> list of (anom, anchor, is_der) for balanced keys
    balanced_t_keys_by_clone: dict[str, dict[tuple[str, ...], list[tuple[str, str, bool]]]] = {}
    # clone -> tkey -> set of breakpoint strings (non-empty)
    t_bp_by_clone: dict[str, dict[tuple[str, ...], set[str]]] = {}

    def t_key_from_str(text: str) -> tuple[str, ...] | None:
        m = re.match(r"^t\(([^)]*)\)", text, re.IGNORECASE)
        if not m:
            return None
        chroms_raw = m.group(1)
        chroms = [c.lstrip("?") for c in chroms_raw.split(";") if c]
        if not chroms:
            return None
        return tuple(sorted(chroms))
    for anom in counts:
        norm = normalize_anomaly(anom)
        base = strip_sign(norm)
        base_core = strip_multiplicity(base)
        if not base_core.startswith(("der", "ider")):
            continue
        t_entries = t_entries_in_str(base_core)
        if not t_entries:
            continue
        if base_core.startswith("ider"):
            ider_t_any[anom] = {k for k, _ in t_entries}
            anchor_match = re.match(r"^ider\(([^)]+)\)", base_core, re.IGNORECASE)
            anchor = anchor_match.group(1) if anchor_match else ""
            is_der = False
        else:
            der_t_any[anom] = {k for k, _ in t_entries}
            anchor_match = re.match(r"^der\(([^)]+)\)", base_core, re.IGNORECASE)
            anchor = anchor_match.group(1) if anchor_match else ""
            is_der = True
        for key, bp in t_entries:
            for clone in clone_map.get(anom, []):
                t_pairs_by_clone.setdefault(clone, {}).setdefault(key, []).append(
                    (anom, anchor, is_der)
                )
                if bp:
                    t_bp_by_clone.setdefault(clone, {}).setdefault(key, set()).add(bp)
        for entry in extract_t_entries(base_core):
            key = entry["key"]
            bp = entry["bp"]
            bucket = der_t_map.setdefault(key, {"bp": bp, "anoms": [], "t_str": entry["t_str"]})
            if bucket["bp"] is None and bp:
                bucket["bp"] = bp
                bucket["t_str"] = entry["t_str"]
            if bp and bucket["bp"] and bp != bucket["bp"]:
                continue
            bucket["anoms"].append(anom)

    for _, bucket in der_t_map.items():
        anoms = bucket["anoms"]
        if len(anoms) < 2:
            continue
        t_str = bucket["t_str"]
        for anom in anoms:
            der_t_by_anom[anom] = t_str

    # équilibrée si au moins un der est présent + un autre der/ider avec un ancrage différent (même clone)
    for clone, tmap in list(t_pairs_by_clone.items()):
        filtered = {}
        for k, items in tmap.items():
            anchors = {anchor for _, anchor, _ in items if anchor}
            has_der = any(is_der for _, _, is_der in items)
            bp_set = t_bp_by_clone.get(clone, {}).get(k, set())
            # équilibrée si 2 ancrages et pas d'incertitude sur les points de cassure
            if len(anchors) >= 2 and has_der:
                bp_has_uncertainty = any("?" in bp for bp in bp_set)
                if bp_has_uncertainty:
                    continue
                # équilibrée si aucun bp explicite ou bp explicite identique
                if len(bp_set) <= 1:
                    filtered[k] = items
        balanced_t_keys_by_clone[clone] = filtered
    explicit_t_keys_by_clone: dict[str, set[tuple[str, ...]]] = {}
    for a in anomalies:
        norm_a = normalize_anomaly(a)
        if not norm_a.startswith("t("):
            continue
        tkey = t_key_from_str(norm_a)
        if not tkey:
            continue
        for clone in clone_map.get(a, []):
            explicit_t_keys_by_clone.setdefault(clone, set()).add(tkey)
    counted_t_keys: dict[str, set[tuple[str, ...]]] = {}
    intertwined_der_scores: dict[str, RuleDecision] = {}

    def anchor_chromosome(anom_str: str) -> str:
        base = strip_multiplicity(strip_sign(normalize_anomaly(anom_str)))
        m = re.match(r"^(?:der|dic)\(([^);]+)", base, re.IGNORECASE)
        return m.group(1).lstrip("?").upper() if m else ""

    def has_intertwined_insertion(anom_str: str, chroms: set[str]) -> bool:
        base = strip_multiplicity(strip_sign(normalize_anomaly(anom_str)))
        return "ins(" in base or len(chroms) >= 3

    def _normalize_breakpoint(bp: str | None) -> str | None:
        if bp is None:
            return None
        return bp.strip().upper()

    def der_breakpoint_signatures(
        anom_str: str, include_uncertain: bool = False
    ) -> set[tuple[str, str]]:
        base = strip_multiplicity(strip_sign(normalize_anomaly(anom_str)))
        signatures: set[tuple[str, str]] = set()

        def add_signature(chrom: str, bp: str | None):
            if not chrom or not bp:
                return
            if not include_uncertain and ("?" in chrom or "?" in bp):
                return
            if "?" not in chrom:
                signatures.add((chrom.upper(), bp.upper()))

        def split_breakpoint_range(text: str | None) -> list[str]:
            bp = _normalize_breakpoint(text)
            if not bp:
                return []
            matches = re.findall(r"[PQ](?:TER|\d+(?:\.\d+)?)", bp)
            return matches or [bp]

        for chroms_raw, bp_raw in re.findall(
            r"t\(([^)]*)\)\(([^)]*)\)", base, re.IGNORECASE
        ):
            chroms = [c.lstrip("?").upper() for c in chroms_raw.split(";") if c]
            bps = [_normalize_breakpoint(bp) for bp in bp_raw.split(";")]
            if len(chroms) != len(bps):
                continue
            for chrom, bp in zip(chroms, bps):
                add_signature(chrom, bp)

        for chroms_raw, bp_raw in re.findall(
            r"ins\(([^)]*)\)\(([^)]*)\)", base, re.IGNORECASE
        ):
            chroms = [c.lstrip("?").upper() for c in chroms_raw.split(";") if c]
            bps = [_normalize_breakpoint(bp) for bp in bp_raw.split(";")]
            if len(chroms) < 2 or len(bps) < 2:
                continue
            add_signature(chroms[0], bps[0])
            for bp in split_breakpoint_range(bps[1]):
                add_signature(chroms[1], bp)

        for expanded in re.findall(r"\(([^()]*(?:->|::)[^()]*)\)", base):
            for chrom, bp in re.findall(
                r"(\d+|X|Y)([pq](?:ter|\d+(?:\.\d+)?))",
                expanded,
                re.IGNORECASE,
            ):
                add_signature(chrom, bp)

        return signatures

    def has_exact_intertwined_breakpoints(
        component_items: list[tuple[str, set[str]]]
    ) -> bool:
        signature_counts: Counter[tuple[str, str]] = Counter()
        for anom, _ in component_items:
            signature_counts.update(der_breakpoint_signatures(anom))
        paired = {sig for sig, count in signature_counts.items() if count >= 2}
        if len(paired) < 3:
            return False
        paired_chroms = {chrom for chrom, _ in paired}
        component_chroms = set().union(*(chroms for _, chroms in component_items))
        return component_chroms.issubset(paired_chroms)

    def has_supported_intertwined_breakpoints(
        component_items: list[tuple[str, set[str]]]
    ) -> bool:
        signature_counts: Counter[tuple[str, str]] = Counter()
        for anom, _ in component_items:
            signature_counts.update(
                der_breakpoint_signatures(anom, include_uncertain=True)
            )
        paired = {sig for sig, count in signature_counts.items() if count >= 2}
        return len(paired) >= 2

    def prepare_intertwined_der_scores():
        by_clone: dict[str, list[tuple[str, set[str]]]] = {}
        for anom in counts:
            base = strip_multiplicity(strip_sign(normalize_anomaly(anom)))
            if not base.startswith(("der", "dic")):
                continue
            chroms = {c for c in get_chromosomes(base) if c != "?"}
            if len(chroms) < 2:
                continue
            for clone in clone_map.get(anom, []):
                by_clone.setdefault(clone, []).append((anom, chroms))

        for _, items in by_clone.items():
            if len(items) < 3:
                continue
            union = set().union(*(chroms for _, chroms in items))
            anchors = {anchor_chromosome(anom) for anom, _ in items}
            if len(union) < 3 or not union.issubset(anchors | union):
                continue
            if len(anchors & union) < 3:
                continue

            remaining = set(range(len(items)))
            while remaining:
                stack = [remaining.pop()]
                component = []
                comp_chroms: set[str] = set()
                while stack:
                    idx = stack.pop()
                    component.append(idx)
                    comp_chroms.update(items[idx][1])
                    linked = [
                        other
                        for other in list(remaining)
                        if items[idx][1].intersection(items[other][1])
                    ]
                    for other in linked:
                        remaining.remove(other)
                        stack.append(other)

                if len(component) < 3 or len(comp_chroms) < 3:
                    continue
                component_items = [items[idx] for idx in component]
                component_anchors = {
                    anchor_chromosome(anom) for anom, _ in component_items
                }
                if len(component_anchors & comp_chroms) < 3:
                    continue
                component_items.sort(
                    key=lambda item: first_index.get(normalize_anomaly(item[0]), float("inf"))
                )
                is_exact = has_exact_intertwined_breakpoints(component_items)
                is_supported = has_supported_intertwined_breakpoints(component_items)
                if not is_exact and not is_supported:
                    continue
                insertion = any(
                    has_intertwined_insertion(anom, chroms)
                    for anom, chroms in component_items
                )
                primary, _ = component_items[0]
                if is_exact:
                    rule_id = (
                        "ISCN.INTERTWINED_DER_BALANCED_INSERTION"
                        if insertion
                        else "ISCN.INTERTWINED_DER_BALANCED_SIMPLE"
                    )
                    score = 2 if insertion else 1
                    explanation = "Dérivés enchevêtrés équilibrés"
                    if insertion:
                        explanation = f"{explanation} (+1 insertion)"
                else:
                    rule_id = "ISCN.INTERTWINED_DER_UNBALANCED"
                    score = 2
                    explanation = (
                        "Dérivés enchevêtrés reliés avec cassures non concordantes"
                    )
                intertwined_der_scores[primary] = RuleDecision(
                    rule_id=rule_id,
                    score=score,
                    explanation=explanation,
                )
                for anom, _ in component_items[1:]:
                    intertwined_der_scores[anom] = RuleDecision(
                        rule_id="ISCN.INTERTWINED_DER_PART",
                        score=0,
                        explanation="Dérivé déjà compté dans le remaniement enchevêtré",
                    )

    prepare_intertwined_der_scores()

    def has_balanced_mirror_der(anom_str: str) -> bool:
        t_keys = der_t_any.get(anom_str, set())
        if not t_keys:
            return False
        chroms = get_chromosomes(strip_multiplicity(strip_sign(normalize_anomaly(anom_str))))
        known_chroms = {chrom for chrom in chroms if chrom != "?"}
        if len(known_chroms) < 2:
            return False
        for clone in clone_map.get(anom_str, []):
            for t_key in t_keys:
                if t_key in balanced_t_keys_by_clone.get(clone, {}):
                    return True
        return False

    def has_der_with_same_t(anom_str: str) -> bool:
        key = t_key_from_str(anom_str)
        if not key:
            return False
        anom_clones = set(clone_map.get(anom_str, []))
        if not anom_clones:
            return False
        for der_anom, der_clones in clone_map.items():
            if "der" not in der_anom:
                continue
            if not anom_clones.intersection(der_clones):
                continue
            base = strip_multiplicity(strip_sign(normalize_anomaly(der_anom)))
            for m in re.finditer(r"t\(([^)]*)\)", base, re.IGNORECASE):
                der_key = tuple(sorted([c.lstrip('?') for c in m.group(1).split(';') if c]))
                if der_key == key:
                    return True
        return False

    def prior_structural_base(anom_str: str) -> str | None:
        norm = normalize_anomaly(anom_str)
        if not norm.startswith("+"):
            return None
        base = strip_multiplicity(strip_sign(norm))
        structural_prefixes = (
            "der",
            "dic",
            "idic",
            "ider",
            "i(",
            "inv",
            "del",
            "dup",
            "ins",
            "add",
            "r(",
        )
        if not base.startswith(structural_prefixes):
            return None
        current_idx = first_index.get(norm, float("inf"))
        for previous in anomalies:
            prev_norm = normalize_anomaly(previous)
            if first_index.get(prev_norm, float("inf")) >= current_idx:
                continue
            prev_base = strip_multiplicity(strip_sign(prev_norm))
            if prev_base == base:
                return base
        return None

    for anom, cnt in counts.items():
        norm = normalize_anomaly(anom)
        base = strip_sign(norm)
        base_core = strip_multiplicity(base)
        cnt_norm = norm_counts[norm]
        eff_cnt = effective_occurrences(anom, cnt, clone_map)

        # a) Constitutionnelles (+Nc) → ISCN = 0
        zeroed_entry = None
        if zeroed_reasons:
            zeroed_entry = zeroed_reasons.get(anom) or zeroed_reasons.get(normalize_anomaly(anom))
        if zeroed_entry:
            override_score, reason = zeroed_entry
            decision = absent_in_repeat_rule_decision("ISCN", override_score, reason)
        else:
            decision = constitutional_rule_decision("ISCN", norm)

        if decision is None:
            decision = apply_rule(
                is_repeat_notation(norm),
                "ISCN.REPEAT_NOTATION",
                0,
                "Anomalies déjà connues dans un autre clone",
            )

        if decision is None and anom in intertwined_der_scores:
            decision = intertwined_der_scores[anom]

        if decision is None and norm.startswith("t("):
            tkey = t_key_from_str(norm)
            clone_keys = set(clone_map.get(anom, []))
            has_pair = False
            is_balanced = False
            if tkey:
                for clone in clone_keys:
                    if tkey in t_pairs_by_clone.get(clone, {}):
                        has_pair = True
                    if tkey in balanced_t_keys_by_clone.get(clone, {}):
                        is_balanced = True
            if tkey and has_pair:
                decision = RuleDecision(
                    rule_id="ISCN.BALANCED_T" if is_balanced else "ISCN.UNBALANCED_T",
                    score=1 if is_balanced else 2,
                    explanation="Translocation équilibrée (dérivés réciproques)" if is_balanced else "Translocation déséquilibrée",
                )
                for clone in clone_keys:
                    counted_t_keys.setdefault(clone, set()).add(tkey)
            elif has_der_with_same_t(anom):
                decision = RuleDecision(
                    rule_id="ISCN.T_ALREADY_IN_DER",
                    score=0,
                    explanation="Translocation déjà comptée dans un dérivé",
                )

        if decision is None:
            if norm in implicit_info:
                info = implicit_info[norm]
                decision = RuleDecision(
                    rule_id="ISCN.IMPLICIT",
                    score=0,
                    explanation=f"{info['reason']} ({info['ref']})",
                )

        if decision is None:
            prior_base = prior_structural_base(anom)
            if prior_base:
                decision = RuleDecision(
                    rule_id="ISCN.STRUCTURAL_GAIN_DUPLICATE",
                    score=1,
                    explanation=(
                        "Duplication d'une anomalie structurale deja decrite "
                        f"({prior_base})"
                    ),
                )

        if decision is None and norm.startswith("+"):
            if base_core.startswith("i(") or base_core.startswith("del"):
                chroms = get_chromosomes(base_core)
                chr_label = chroms and sorted(chroms)[0] or ""
                if base_core.startswith("i("):
                    expl = "Équivalence sémantique: +i(...) = +chr + i(...)"
                elif base_core.startswith("del"):
                    expl = "Équivalence sémantique: +del(...) = +chr + del(...)"
                if chr_label:
                    expl = f"{expl} (chr {chr_label})"
                decision = RuleDecision(
                    rule_id="ISCN.SEMANTIC_PLUS_ISO" if base_core.startswith("i(") else "ISCN.SEMANTIC_PLUS_DEL",
                    score=2,
                    explanation=expl,
                )

        if decision is None:
            decision = apply_rule(
                is_marker_anomaly(norm),
                "ISCN.MAR",
                1,
                "Marqueur",
            )

        if decision is None:
            decision = apply_rule(
                base_core.startswith('dic'),
                "ISCN.DICENTRIC",
                2,
                "Chromosome dicentrique",
            )

        if decision is None and base_core.startswith('der'):
            t_keys = der_t_any.get(anom, set())
            clone_keys = clone_map.get(anom, [])
            chosen_key = None
            chosen_clone = None
            is_balanced = False
            for clone in clone_keys:
                for k in t_keys:
                    if k in t_pairs_by_clone.get(clone, {}):
                        chosen_key = k
                        chosen_clone = clone
                        is_balanced = k in balanced_t_keys_by_clone.get(clone, {})
                        break
                if chosen_key:
                    break
            if chosen_key and chosen_clone:
                # Si second chromosome incertain, ne pas surcoter via la t
                chroms = get_chromosomes(base_core)
                known_chroms = {chrom for chrom in chroms if chrom != "?"}
                if len(chosen_key) < 2 or len(known_chroms) < 2 or "?" in chroms:
                    decision = None
                elif chosen_key in explicit_t_keys_by_clone.get(chosen_clone, set()):
                    decision = RuleDecision(
                        rule_id="ISCN.DER_T_COUNTED",
                        score=0,
                        explanation="Dérivé associé à une translocation comptée",
                    )
                elif chosen_key not in counted_t_keys.get(chosen_clone, set()):
                    decision = RuleDecision(
                        rule_id="ISCN.T_VIA_DER_BALANCED" if is_balanced else "ISCN.T_VIA_DER_UNBALANCED",
                        score=1 if is_balanced else 2,
                        explanation="Translocation équilibrée (comptée via dérivé)" if is_balanced else "Translocation déséquilibrée (comptée via dérivé)",
                    )
                    counted_t_keys.setdefault(chosen_clone, set()).add(chosen_key)
                else:
                    decision = RuleDecision(
                        rule_id="ISCN.DER_T_COUNTED",
                        score=0,
                        explanation="Dérivé associé à une translocation comptée",
                    )

        if decision is None and base_core.startswith('der'):
            if anom in der_t_by_anom and not has_extra_der_components(base_core):
                decision = RuleDecision(
                    rule_id="ISCN.DER_BALANCED_T",
                    score=0,
                    explanation=f"Dérivé issu d'une translocation équilibrée ({der_t_by_anom[anom]})",
                )

        if decision is None and base_core.startswith('der'):
            chroms = get_chromosomes(base_core)
            known = count_known_chromosomes(chroms)
            has_unknown = '?' in chroms
            if is_der_without_breakpoints(base_core):
                decision = RuleDecision(
                    rule_id="ISCN.DER_NO_BREAKPOINT",
                    score=1,
                    explanation="Chromosome dérivé sans point de cassure détaillé",
                )
            elif known >= 2:
                if has_unknown and '?' in base_core:
                    decision = RuleDecision(
                        rule_id="ISCN.DER_MULTI_UNCERTAIN",
                        score=2,
                        explanation="Chromosome dérivé impliquant plusieurs chromosomes avec des imprécisions",
                    )
                else:
                    decision = RuleDecision(
                        rule_id="ISCN.DER_MULTI",
                        score=2,
                        explanation="Chromosome dérivé impliquant plusieurs chromosomes",
                    )
            elif has_unknown:
                decision = RuleDecision(
                    rule_id="ISCN.DER_UNCERTAIN_SECOND",
                    score=1,
                    explanation="Chromosome dérivé sans certitude sur un second chromosome",
                )
            else:
                decision = RuleDecision(
                    rule_id="ISCN.DER_SAME_CHR",
                    score=1,
                    explanation="Chromosome dérivé issu du même chromosome",
                )

        # Si un der contient une anomalie structurale additionnelle, on ajoute +1.
        # Liste restreinte (à élargir si besoin) : add, inv, ins, del, dup.
        # Si un der est préfixé par '+', on ajoute +1 pour le gain du chromosome d'ancrage.
        extra_component = None
        extra_der_gain = False
        if base_core.startswith("der"):
            extra_match = re.search(r"(add|inv|ins|del|dup)\(", base_core)
            if extra_match and has_balanced_mirror_der(anom):
                extra_component = extra_match.group(1)
            if (
                norm.startswith("+")
                and (
                    decision is None
                    or decision.rule_id != "ISCN.STRUCTURAL_GAIN_DUPLICATE"
                )
            ):
                extra_der_gain = True

        if decision is None and norm.startswith(("+", "-")):
            decision = single_chr_deseq_rule_decision(norm, cnt_norm)
            if decision is None:
                decision = apply_rule(
                    is_complex_multichr_deseq(base_core),
                    "ISCN.COMPLEX_MULTI_CHR",
                    2,
                    "Déséquilibre multichromosomique complexe",
                )
            if decision is None:
                decision = apply_rule(
                    is_unbalanced_translocation(base_core),
                    "ISCN.UNBALANCED_TRANSLOCATION",
                    2,
                    "Translocation déséquilibrée",
                )
            if decision is None:
                decision = RuleDecision(
                    rule_id="ISCN.GAIN_SIMPLE" if norm.startswith("+") else "ISCN.LOSS_SIMPLE",
                    score=1,
                    explanation="-",
                )

        if decision is None:
            decision = single_chr_deseq_rule_decision(norm, cnt_norm)
            if decision is None:
                decision = apply_rule(
                    is_complex_multichr_deseq(base_core),
                    "ISCN.COMPLEX_MULTI_CHR",
                    2,
                    "Déséquilibre multichromosomique complexe",
                )
            if decision is None:
                decision = apply_rule(
                    is_unbalanced_translocation(base_core),
                    "ISCN.UNBALANCED_TRANSLOCATION",
                    2,
                    "Translocation déséquilibrée",
                )
            if decision is None:
                decision = RuleDecision(
                    rule_id="ISCN.OTHER_STANDARD",
                    score=1,
                    explanation="-",
                )

        score = decision.score
        explication = decision.explanation
        score_multiplier = eff_cnt
        if decision.rule_id.startswith("ISCN.ABSENT_IN_REPEAT") and decision.score == 1:
            score_multiplier = 1
        if (
            decision.rule_id == "ISCN.SINGLE_CHR_GAIN_REPEAT"
            and norm.startswith("+")
            and cnt_norm > 1
        ):
            score_multiplier = 1
            explication = f"{explication} (gain répété : {anom})"
        score = decision.score * score_multiplier
        if (
            decision.rule_id.startswith("ISCN.ABSENT_IN_REPEAT")
            and "Perte chromosomique implicite" in explication
        ):
            score = 1
        if extra_component:
            score += 1
            explication = f"{explication} (+1 pour {extra_component})"
        if extra_der_gain:
            score += 1
            explication = f"{explication} (+1 pour +chr)"

        total += score
        explication = append_uncertainty_note(anom, explication)

        rows.append({
            "Anomalie": anom,
            "Type": type_anomalie(norm),
            "Explication": explication,
            "Occurrences": eff_cnt,
            "Clones": ", ".join(clone_map.get(anom, [])),
            "Score ISCN 2024": score,
        })
        rule_ids[anom] = decision.rule_id
        rule_expls[anom] = decision.explanation

    # Ligne de totaux
    rows.append({
        "Anomalie": "TOTAL",
        "Type": "",
        "Explication": "",
        "Occurrences": "",
        "Clones": "",
        "Score ISCN 2024": total,
    })

    return pd.DataFrame(rows), total, rule_ids, rule_expls

# Fonction pour analyser une formule caryotypique
# =========================
# Orchestration
# =========================
def deduplicate_inter_clones(entries):
    """Déduplication inter-clones sans modifier la logique existante."""

    first_clone_by_norm: dict[str, str] = {}
    scorable_entries: list[dict[str, str]] = []
    clone_details_map: dict[str, dict[str, dict[str, object]]] = {}
    marker_suppressed: dict[str, set[str]] = {}
    marker_ref_counts: dict[str, int] = {}
    zeroed_reasons: dict[str, tuple[int, str]] = {}

    def unit_entry(entry: dict) -> dict:
        """Retourne une occurrence unitaire pour éviter une double multiplicité."""

        copied = dict(entry)
        copied["count"] = 1
        return copied

    for entry in entries:
        raw_anom = entry["anomaly"]
        norm = normalize_anomaly(raw_anom)
        clone = entry["clone"]
        count = int(entry.get("count", 1))
        is_repeat_clone = bool(entry.get("repeat_clone"))
        is_marker = is_marker_anomaly(norm)
        key_for_ref = strip_multiplicity(norm) if is_marker else raw_anom

        if key_for_ref not in first_clone_by_norm:
            first_clone_by_norm[key_for_ref] = clone
        ref_clone = first_clone_by_norm[key_for_ref]
        is_reference = clone == ref_clone

        clone_bucket = clone_details_map.setdefault(entry["anomaly"], {})

        def get_bucket(key, is_ref, reason):
            bucket_entry = clone_bucket.setdefault(key, {
                "label": display_clone_label(clone),
                "is_reference": is_ref,
                "reason": reason,
                "count": 0,
                "score_override": None,
            })
            return bucket_entry

        suppressed_one = False
        scorable_count = count if is_reference else 0
        if is_repeat_clone and is_negative_anomaly(raw_anom):
            stripped = strip_negative(raw_anom)
            stripped_norm = normalize_anomaly(strip_sign(stripped))

            def _matches_previous() -> bool:
                if is_marker_anomaly(stripped_norm):
                    target = strip_sign(stripped_norm)
                    for prev in first_clone_by_norm.keys():
                        if is_marker_anomaly(prev) and strip_sign(prev) == target:
                            return True
                    return False
                if stripped_norm.startswith(("add", "del", "dup", "ins", "inv")):
                    try:
                        base = stripped_norm.split("(")[0] + "(" + stripped_norm.split("(")[1].split(")")[0] + ")"
                    except IndexError:
                        return False
                    for prev in first_clone_by_norm.keys():
                        if prev.startswith(base):
                            return True
                    return False
                if stripped_norm.startswith(("der", "dic")):
                    m = re.match(r"^(der|dic)\(([^)]*)\)(?:t\(([^)]*)\))?", stripped_norm, re.IGNORECASE)
                    if not m:
                        return False
                    kind, der_chr, t_chr = m.groups()
                    for prev in first_clone_by_norm.keys():
                        pm = re.match(rf"^{kind}\(([^)]*)\)(?:t\(([^)]*)\))?", prev, re.IGNORECASE)
                        if not pm:
                            continue
                        p_der, p_t = pm.groups()
                        if p_der == der_chr and (t_chr is None or p_t == t_chr):
                            return True
                    return False
                return stripped_norm in first_clone_by_norm

            if _matches_previous():
                reason = f"Anomalie indiquée absente dans {display_clone_label(clone)} (déjà comptée)"
                score_override = 0
                if stripped_norm.startswith(("der", "dic", "add")):
                    reason = f"Perte chromosomique implicite dans {display_clone_label(clone)}"
                    score_override = 1
                elif stripped_norm.startswith("r("):
                    if re.search(r"r\(\s*\d+", stripped_norm):
                        reason = f"Perte chromosomique implicite dans {display_clone_label(clone)}"
                        score_override = 1
                bucket_entry = get_bucket(clone, False, reason)
                bucket_entry["count"] = bucket_entry.get("count", 0) + count
                bucket_entry["score_override"] = score_override
                zeroed_reasons[raw_anom] = (score_override, reason)
                zeroed_reasons[normalize_anomaly(raw_anom)] = (score_override, reason)
                for _ in range(count):
                    scorable_entries.append(unit_entry(entry))
                continue
        if is_reference and is_marker:
            marker_ref_counts.setdefault(key_for_ref, count)
        if not is_reference and is_marker:
            ref_count = marker_ref_counts.get(key_for_ref, 0)
            suppressed = marker_suppressed.setdefault(key_for_ref, set())
            if clone not in suppressed and ref_count > 0:
                suppressed.add(clone)
                suppressed_one = True
            scorable_count = max(count - ref_count, 0)

        if scorable_count > 0:
            reason = ""
            if suppressed_one:
                ref_count = marker_ref_counts.get(key_for_ref, 1)
                reason = f"{ref_count} marqueur déjà compté dans {display_clone_label(ref_clone)}"
            bucket_entry = get_bucket(clone, True, reason)
            display_increment = 1 if is_marker else scorable_count
            bucket_entry["count"] = bucket_entry.get("count", 0) + display_increment
            if not is_reference and is_marker:
                bucket_entry["is_reference"] = True
            for _ in range(scorable_count):
                scorable_entries.append(unit_entry(entry))

    clone_details_ready = {
        anom: list(clone_dict.values()) for anom, clone_dict in clone_details_map.items()
    }

    return scorable_entries, clone_details_ready, zeroed_reasons


def analyser_formule(formule, debug: bool = False):
    """
    Analyse une formule caryotypique et retourne:
    - Le DataFrame des anomalies détectées
    - Un dictionnaire de scores totaux (ISCN et Jondreville)
    - Une erreur éventuelle
    """
    try:
        expanded_formula, condensed_changed = expand_condensed_formula(formule)
        anomalies, clone_map, entries = parse_caryotype(expanded_formula)

        # Déduplication inter-clones: seule la première apparition d'une
        # anomalie (normalisée) est comptabilisée, les autres clones obtiennent
        # un score nul avec justification.
        scorable_entries, clone_details_ready, zeroed_reasons = deduplicate_inter_clones(entries)
        scorable_anomalies = [entry["anomaly"] for entry in scorable_entries]
        scoring_clone_map = {}
        for entry in scorable_entries:
            scoring_clone_map.setdefault(entry["anomaly"], []).append(entry["clone"])

        df_iscn, total_iscn, rule_id_iscn, rule_expl_iscn = calcul_score_iscn(
            scorable_anomalies, clone_map, zeroed_reasons, scoring_clone_map
        )
        jondroville_scores, jondroville_explanations, total_jondroville, rule_id_jon, rule_expl_jon = calcul_score_jondroville(
            scorable_anomalies, clone_map, scorable_entries, zeroed_reasons
        )

        df_iscn["Score Jondreville 2020"] = df_iscn["Anomalie"].apply(
            lambda anom: total_jondroville if anom == "TOTAL" else jondroville_scores.get(anom, 0)
        )

        df_iscn["Explication Jondreville 2020"] = df_iscn["Anomalie"].apply(
            lambda anom: "" if anom == "TOTAL" else jondroville_explanations.get(anom, "Anomalie non constitutionnelle")
        )

        df_iscn["CloneDetails"] = df_iscn["Anomalie"].apply(
            lambda anom: [] if anom == "TOTAL" else clone_details_ready.get(anom, [])
        )

        if debug:
            df_iscn["RuleID_ISCN"] = df_iscn["Anomalie"].apply(
                lambda anom: "" if anom == "TOTAL" else rule_id_iscn.get(anom, "")
            )
            df_iscn["RuleID_Jon"] = df_iscn["Anomalie"].apply(
                lambda anom: "" if anom == "TOTAL" else rule_id_jon.get(anom, "")
            )
            df_iscn["RuleExplanation_ISCN"] = df_iscn["Anomalie"].apply(
                lambda anom: "" if anom == "TOTAL" else rule_expl_iscn.get(anom, "")
            )
            df_iscn["RuleExplanation_Jon"] = df_iscn["Anomalie"].apply(
                lambda anom: "" if anom == "TOTAL" else rule_expl_jon.get(anom, "")
            )

        return df_iscn, {
            "iscn": total_iscn,
            "jondroville": total_jondroville,
            "formule_originale": formule if condensed_changed else "",
            "formule_equivalente": expanded_formula if condensed_changed else "",
        }, None
    except Exception as e:
        return None, {"iscn": 0, "jondroville": 0}, f"Erreur lors de l'analyse de la formule: {str(e)}"
