import re
from collections import Counter
import pandas as pd
from dataclasses import dataclass

# =========================
# Parsing
# =========================

# Extraction des numéros de chromosome dans une anomalie ISCN
def get_chromosomes(anom):
    """Retourne l'ensemble des chromosomes impliqués dans ``anom``.

    La fonction détecte les numéros apparaissant :
    - juste après les mots clés (der, del, dup, t, ...)
    - dans la seconde parenthèse des notations ``der(...)`` (apès les flèches)
    - précédés d'un ``?`` comme dans ``t(?1;17)``
    """

    nums: set[str] = set()

    # 1) Numéros directement après der(...), t(...), etc.
    for m in re.finditer(r'(?:der|dic|del|dup|ins|t|i|ider|idic|r)\((\??[0-9;?]+)', anom):
        raw = m.group(1)
        # Se limiter à la partie numérique avant un ")" ou une nouvelle parenthèse
        raw = re.split(r'[)()]', raw)[0]
        for num in raw.split(';'):
            cleaned = num.lstrip('?')
            if cleaned:
                nums.add(cleaned)
            elif '?' in num:
                nums.add('?')

    # 2) Numéros mentionnés dans la seconde parenthèse des der(...)
    for _, second in re.findall(r'der\(([^)]*)\)\(([^)]*)\)', anom):
        for n in re.findall(r'\??(\d+)(?=[pq])', second):
            nums.add(n.lstrip('?'))
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
    """Nombre d'occurrences à scorer après contexte de ploidie."""

    norm = normalize_anomaly(anom)
    multiplicity = anomaly_occurrences(norm)
    if (
        multiplicity == 2
        and is_tetraploid_context(anom, clone_map)
        and tetraploid_clone_uses_global_multiplicity(anom, clone_map)
    ):
        return max(count // multiplicity, 1)
    return count


def is_tetraploid_context(anom: str, clone_map: dict[str, list[str]]) -> bool:
    """Indique si l'anomalie appartient à un clone en tétraploïdie."""

    tetrap_clones = set(clone_map.get("Tetraploidy", []))
    if not tetrap_clones:
        return False
    anom_clones = set(clone_map.get(anom, []))
    return bool(tetrap_clones & anom_clones)


def tetraploid_clone_uses_global_multiplicity(
    anom: str, clone_map: dict[str, list[str]]
) -> bool:
    """Indique si le x2 est une notation globale du clone tetraploide.

    Si le clone contient aussi des anomalies non suffixees x2, le x2 restant est
    interprete comme une vraie multiplicite locale de cette anomalie.
    """

    tetrap_clones = set(clone_map.get("Tetraploidy", []))
    anom_clones = set(clone_map.get(anom, [])) & tetrap_clones
    if not anom_clones:
        return False

    for clone in anom_clones:
        clone_anomalies = [
            other
            for other, clones in clone_map.items()
            if clone in clones and other not in {"Tetraploidy", "Triploidy"}
        ]
        if clone_anomalies and all(
            anomaly_occurrences(normalize_anomaly(other)) == 2
            for other in clone_anomalies
        ):
            return True
    return False


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
                # Chromosomes juste apres der(...)
                chrs = extract_chr_ids(m.group(1))
                # Ajouter egalement les partenaires de la/les translocations t(...)
                for t in re.finditer(r"t\(([0-9?;]+)\)", an):
                    chrs.update(extract_chr_ids(t.group(1)))
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
class RuleDecision:
    rule_id: str
    score: int
    explanation: str


def apply_rule(condition: bool, rule_id: str, score: int, explanation: str) -> RuleDecision | None:
    if condition:
        return RuleDecision(rule_id=rule_id, score=score, explanation=explanation)
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
            filtered.append(anom)
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
        is_constitutional, const_expl = constitutional_status(norm)
        zeroed_entry = None
        if zeroed_reasons:
            zeroed_entry = zeroed_reasons.get(anom) or zeroed_reasons.get(normalize_anomaly(anom))
        if zeroed_entry:
            override_score, reason = zeroed_entry
            decision = RuleDecision(
                rule_id="JON.ABSENT_IN_REPEAT",
                score=override_score,
                explanation=reason,
            )
        else:
            decision = apply_rule(
                is_constitutional,
                "JON.CONSTITUTIONAL",
                0,
                const_expl,
            )
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
        if decision.rule_id == "JON.ABSENT_IN_REPEAT":
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
):
    """Calcule le détail des scores selon la grille ISCN 2024."""

    counts = Counter(anomalies)
    norm_counts = Counter(normalize_anomaly(a) for a in anomalies)
    implicit_info = detect_implicit_anomalies(anomalies, clone_map)
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
        is_constitutional, const_expl = constitutional_status(norm)
        zeroed_entry = None
        if zeroed_reasons:
            zeroed_entry = zeroed_reasons.get(anom) or zeroed_reasons.get(normalize_anomaly(anom))
        if zeroed_entry:
            override_score, reason = zeroed_entry
            decision = RuleDecision(
                rule_id="ISCN.ABSENT_IN_REPEAT",
                score=override_score,
                explanation=reason,
            )
        else:
            decision = apply_rule(
                is_constitutional,
                "ISCN.CONSTITUTIONAL",
                0,
                const_expl,
            )

        if decision is None:
            decision = apply_rule(
                is_repeat_notation(norm),
                "ISCN.REPEAT_NOTATION",
                0,
                "Anomalies déjà connues dans un autre clone",
            )

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
                    rule_id="ISCN.SEMANTIC_PLUS_STRUCT",
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
                if len(chosen_key) < 2 or "?" in chroms:
                    decision = None
                elif chosen_key in explicit_t_keys_by_clone.get(chosen_clone, set()):
                    decision = RuleDecision(
                        rule_id="ISCN.DER_T_COUNTED",
                        score=0,
                        explanation="Dérivé associé à une translocation comptée",
                    )
                elif chosen_key not in counted_t_keys.get(chosen_clone, set()):
                    decision = RuleDecision(
                        rule_id="ISCN.T_VIA_DER",
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
            decision = apply_rule(
                is_single_chr_deseq(norm, cnt_norm),
                "ISCN.SINGLE_CHR_DESEQ",
                2,
                "Déséquilibre unichromosomique",
            )
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
                    rule_id="ISCN.GAIN_LOSS_SIMPLE",
                    score=1,
                    explanation="-",
                )

        if decision is None:
            decision = apply_rule(
                is_single_chr_deseq(norm, cnt_norm),
                "ISCN.SINGLE_CHR_DESEQ",
                2,
                "Déséquilibre unichromosomique",
            )
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
        if decision.rule_id == "ISCN.ABSENT_IN_REPEAT" and decision.score == 1:
            score_multiplier = 1
        if (
            decision.rule_id == "ISCN.SINGLE_CHR_DESEQ"
            and norm.startswith("+")
            and cnt_norm > 1
        ):
            score_multiplier = 1
            explication = f"{explication} (gain répété : {anom})"
        score = decision.score * score_multiplier
        if (
            decision.rule_id == "ISCN.ABSENT_IN_REPEAT"
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
                    scorable_entries.append(entry)
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
                scorable_entries.append(entry)

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

        df_iscn, total_iscn, rule_id_iscn, rule_expl_iscn = calcul_score_iscn(
            scorable_anomalies, clone_map, zeroed_reasons
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
