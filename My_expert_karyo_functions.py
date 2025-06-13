import re
from collections import Counter
import pandas as pd

# Extraction des numéros de chromosome dans une anomalie ISCN
def get_chromosomes(anom):
    """
    Retourne l'ensemble des chromosomes impliqués dans l'anomalie,
    basé sur les notations avant chaque parenthèse de type der, del, dup, ins, t, i, ider, idic.
    """
    nums = set()
    for m in re.finditer(r'(?:der|dic|del|dup|ins|t|i|ider|idic|r)\(([0-9;]+)', anom):
        for num in m.group(1).split(';'):
            nums.add(num)
    return nums

# Parsing de la formule karyotypique
def parse_caryotype(chaine_iscn):
    """
    Parse une chaîne ISCN avec clones séparés par '/'.
    Renvoie la liste plate des anomalies et un dict {anom: [clones]}.
    Gère aussi les anomalies de ploidie (≠46).
    """
    # Remove all whitespace for robust parsing
    chaine_iscn = re.sub(r"\s+", "", chaine_iscn)
    
    anomalies = []
    clone_map = {}
    clones = [re.sub(r"\[.*?\]", "", c) for c in chaine_iscn.split('/')]
    for idx, clone in enumerate(clones, start=1):
        parts = [p.strip().strip('.') for p in clone.split(',') if p.strip()]
        # Détection de la ploidie
        try:
            total = int(re.sub(r"\D", "", parts[0]))
            if total != 46:
                if total == 92:
                    pl = 'Tetraploidy'
                elif total == 69:
                    pl = 'Triploidy'

                anomalies.append(pl)
                clone_map.setdefault(pl, []).append(f"clone{idx}")
        except Exception:
            pass
        # Extraction des anomalies structurelles
        for an in parts[2:]:
            anomalies.append(an)
            clone_map.setdefault(an, []).append(f"clone{idx}")
    return anomalies, clone_map

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
    pattern = r'^t\(\d+(?:;\d+)+\)\(.+\)$'
    return bool(re.match(pattern, anom)) and 'der' not in anom and '+' not in anom and '-' not in anom

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
    pattern = r'^ins\(\d+(?:;\d+)+\)\(.+\)$'
    return bool(re.match(pattern, anom)) and 'der' not in anom and '+' not in anom and '-' not in anom

# Détection des anomalies multichromosomiques déséquilibrées pour 2 points
def is_complex_multichr_deseq(anom):
    """
    Détecte les anomalies multichromosomiques déséquilibrées (≥2 chromosomes) pour 2 points.
    Les chromosomes dérivés ("der") sont considérés complexes par définition
    même si un seul numéro est explicitement indiqué.
    Renvoie False si un seul chromosome impliqué.
    """
    # Les "der" sont toujours multichromosomiques
    if anom.startswith('der'):
        return True

    chroms = get_chromosomes(anom)
    # si un seul chromosome impliqué -> pas multi-chromosomique déséquilibrée
    if len(chroms) <= 1:
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
    if is_complex_multichr_deseq(anom):
        return 'Multichromosomique déséquilibrée'
    if is_balanced_translocation(anom):
        return 'Translocation équilibrée'
    if is_unbalanced_translocation(anom):
        return 'Translocation déséquilibrée'
    if is_balanced_insertion(anom):
        return 'Insertion équilibrée'
    if anom == '<2n>':
        return 'Ploidy'
    if '~' in anom:
        return 'Pléiade chromosomique'
    if anom == '+mar':
        return 'Chromosome marqueur'
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
    if anom.startswith('t('):
        return 'Translocation'
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
    return 'Autre'

# Calcul des scores
def normalize_anomaly(anom: str) -> str:
    """Normalise une anomalie pour le scoring.

    - Supprime un éventuel point d'interrogation en début d'anomalie
      ("?dic" -> "dic").
    """
    norm = anom.lstrip('?')
    return norm


def detect_implicit_anomalies(anomalies):
    """Détecte les anomalies implicites et renvoie un dict.

    Le dict a pour clé l'anomalie normalisée et pour valeur un
    dictionnaire avec la clef ``reason`` décrivant la cause et ``ref``
    l'anomalie de référence à afficher entre parenthèses.
    """
    norm_counts = Counter(normalize_anomaly(a) for a in anomalies)
    # mappage normalisé -> version originale pour l'affichage
    norm_to_orig = {}
    for a in anomalies:
        norm = normalize_anomaly(a)
        norm_to_orig.setdefault(norm, a)

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
    multi_der = {}
    for an in norm_counts:
        if an.startswith(('der', 'dic')):
            m = re.match(r"^(?:der|dic)\(([0-9;]+)\)", an)
            if m:
                # Chromosomes juste apres der(...)
                chrs = set(m.group(1).split(';'))
                # Ajouter egalement les partenaires de la/les translocations t(...)
                for t in re.finditer(r"t\(([0-9;]+)\)", an):
                    chrs.update(t.group(1).split(';'))
                if len(chrs) > 1:
                    for c in chrs:
                        multi_der.setdefault(c, []).append(an)

    for an in norm_counts:
        if an.startswith(('+', '-')):
            num = re.sub(r"\D", "", an)
            if num in multi_der:
                ref_norm = multi_der[num][0]
                ref = norm_to_orig.get(ref_norm, ref_norm)
                implicit[an] = {"reason": "Gain/perte implicite", "ref": ref}

    # 3) Répétitions de la même anomalie (même chromosome)
    base_pattern = re.compile(
        r'^(?:der|dic|del|dup|ins|t|i|ider|idic|r|add)\([0-9;]+\)'
    )
    base_map = {}
    for a in anomalies:
        norm = normalize_anomaly(a)
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

    return implicit


def calcul_scores(anomalies, clone_map):
    """
    Calcule les scores selon deux méthodes:
    
    Jondreville 2020 : 1 point par anomalie (toujours)
    ISCN 2024     :
      - 0 pt pour anomalies implicites (dérivés implicites & constitutionnelles)
      - 2 pts pour déséquilibres unichr/multichr ou translocations déséquilibrées
      - 1 pt pour anomalies standard
    """
    counts = Counter(anomalies)
    norm_counts = Counter(normalize_anomaly(a) for a in anomalies)

    implicit_info = detect_implicit_anomalies(anomalies)

    rows = []
    total_j = total_i = 0

    for anom, cnt in counts.items():
        score_j = 1  # Jondreville = 1 pour toutes
        norm = normalize_anomaly(anom)
        cnt_norm = norm_counts[norm]

        # a) Constitutionnelles (+Nc) → ISCN = 0
        if re.match(r"^\+\d+c$", norm):
            score_i = 0
            explication = "Anomalie constitutionnelle (0 point)"

        # b) Anomalies détectées comme implicites
        elif norm in implicit_info:
            info = implicit_info[norm]
            score_i = 0
            explication = f"{info['reason']} ({info['ref']}) (0 point)"

        # c) Gains/pertes simples (analyse standard si non implicite)
        elif norm.startswith(("+", "-")):
            if is_single_chr_deseq(norm, cnt_norm):
                score_i = 2
                explication = "Déséquilibre unichromosomique (2 points)"
            elif is_complex_multichr_deseq(norm):
                score_i = 2
                explication = "Déséquilibre multichromosomique complexe (2 points)"
            elif is_unbalanced_translocation(norm):
                score_i = 2
                explication = "Translocation déséquilibrée (2 points)"
            else:
                score_i = 1
                explication = "Anomalie standard (1 point)"

        # d) Chromosomes dicentriques → 2 points
        elif norm.startswith('dic'):
            score_i = 2
            explication = "Chromosome dicentrique (2 points)"

        # e) Toutes les autres anomalies → scoring standard
        else:
            if is_single_chr_deseq(norm, cnt_norm):
                score_i = 2
                explication = "Déséquilibre unichromosomique (2 points)"
            elif is_complex_multichr_deseq(norm):
                score_i = 2
                explication = "Déséquilibre multichromosomique complexe (2 points)"
            elif is_unbalanced_translocation(norm):
                score_i = 2
                explication = "Translocation déséquilibrée (2 points)"
            else:
                score_i = 1
                explication = "Anomalie standard (1 point)"

        total_j += score_j
        total_i += score_i

        rows.append({
            "Anomalie": anom,
            "Type": type_anomalie(norm),
            "Explication": explication,
            "Occurrences": cnt,
            "Clones": ", ".join(clone_map.get(anom, [])),
            "Score Jondreville 2020": score_j,
            "Score ISCN 2024": score_i,
        })

    # Ligne de totaux
    rows.append({
        "Anomalie": "TOTAL",
        "Type": "",
        "Explication": "",
        "Occurrences": "",
        "Clones": "",
        "Score Jondreville 2020": total_j,
        "Score ISCN 2024": total_i,
    })

    return pd.DataFrame(rows), total_i

# Fonction pour analyser une formule caryotypique
def analyser_formule(formule):
    """
    Analyse une formule caryotypique et retourne:
    - Le DataFrame des anomalies détectées
    - Le score total
    - Une erreur éventuelle
    """
    try:
        anomalies, clone_map = parse_caryotype(formule)
        df, total = calcul_scores(anomalies, clone_map)
        return df, total, None
    except Exception as e:
        return None, 0, f"Erreur lors de l'analyse de la formule: {str(e)}"