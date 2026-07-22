import streamlit as st
import pandas as pd
import re
import base64
import io
import html
import math
import uuid
import openpyxl
import streamlit.components.v1 as components
from pathlib import Path
from My_expert_karyo_functions import (
    analyser_formule,
    get_canonical_rule_catalog_dataframe,
    get_rule_catalog_dataframe,
    get_rule_path,
)

LOCAL_TEST_FILENAME = "comptage_local_MYC.xlsx"
LOCAL_TEST_PATH = Path(__file__).resolve().parent / LOCAL_TEST_FILENAME
TEST_SHEET_URL = (
    "https://docs.google.com/spreadsheets/d/"
    "1_QPAzEr3PNHaNVu8Qxuv_hWCUoBQqqRcNSnX-Q0Egm8/edit?gid=98233212#gid=98233212"
)
RULE_CATALOG_SHEET_URL = "http://docs.google.com/spreadsheets/d/1MkwGWtuRU53fuaZ61RejAapUZ1NIeZcIJwMuv7OaY4w/edit?gid=292168720#gid=292168720"
RULE_CATALOG_REFERENCE_PATH = Path(__file__).resolve().parent / "rules_catalog_reference.csv"

# Configuration de la page
st.set_page_config(
    page_title="Analyseur de Caryotypes",
    page_icon="🧬",
    layout="wide"
)

# Utilitaire pour charger une feuille Google Sheets publique
def load_google_sheet(url_or_id):
    """Charge une feuille Google Sheets publique en DataFrame.

    Parameters
    ----------
    url_or_id : str
        L'URL complète ou l'identifiant de la feuille.

    Returns
    -------
    tuple[pd.DataFrame | None, str | None]
        Le DataFrame chargé, ou ``None`` en cas d'erreur, ainsi qu'un message
        d'erreur éventuel.
    """
    if not url_or_id:
        return None, "URL vide"

    # Extraire l'identifiant du document
    match = re.search(r"/d/([A-Za-z0-9_-]+)", url_or_id)
    sheet_id = match.group(1) if match else url_or_id

    # Extraire l'identifiant de feuille (gid) si présent
    gid_match = re.search(r"gid=([0-9]+)", url_or_id)
    gid = gid_match.group(1) if gid_match else "0"

    csv_url = (
        f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq"
        f"?tqx=out:csv&gid={gid}"
    )

    try:
        df = pd.read_csv(csv_url)
        return df, None
    except Exception as e:
        return None, str(e)


def detect_input_columns(df):
    formule_col = None
    for col in df.columns:
        if str(col).strip().lower() == 'formule':
            formule_col = col
            break

    count_i_col = None
    count_j_col = None
    for col in df.columns:
        norm = str(col).strip().lower()
        if norm in ('count_i', 'count_iscn'):
            count_i_col = col
        elif norm in ('count_j', 'count_jon'):
            count_j_col = col

    return formule_col, count_i_col, count_j_col


def normalize_reference_value(value):
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


def floor_percent(ok: int, total: int) -> float | None:
    if not total:
        return None
    return math.floor((ok / total * 100) * 100) / 100


def format_percent(value: float) -> str:
    return f"{value:.2f}".replace(".", ",")


def summarize_matches(match_details):
    """Calcule les mêmes compteurs de concordance pour l'aperçu et le bilan."""

    summary = {}
    for system in ("iscn", "jon"):
        total = sum(1 for match in match_details if match.get(system) is not None)
        ok = sum(1 for match in match_details if match.get(system) is True)
        summary[system] = {
            "ok": ok,
            "total": total,
            "percent": floor_percent(ok, total),
        }
    return summary


def compute_match_preview(df):
    formule_col, count_i_col, count_j_col = detect_input_columns(df)
    if formule_col is None:
        return None

    errors = 0
    match_details = []

    for _, row in df.iterrows():
        _, totals, error = analyser_formule(row[formule_col])
        match_detail = {"iscn": None, "jon": None}
        if error:
            errors += 1
            match_details.append(match_detail)
            continue

        if count_i_col:
            ref_i = normalize_reference_value(row[count_i_col])
            if ref_i is not None:
                match_detail["iscn"] = totals.get("iscn") == ref_i
        if count_j_col:
            ref_j = normalize_reference_value(row[count_j_col])
            if ref_j is not None:
                match_detail["jon"] = totals.get("jondroville") == ref_j
        match_details.append(match_detail)

    summary = summarize_matches(match_details)
    summary["errors"] = errors
    return summary


@st.cache_data(show_spinner=False, ttl=300)
def load_google_sheet_preview_v2(url):
    df, err = load_google_sheet(url)
    if err:
        return None, err
    return compute_match_preview(df), None


@st.cache_data(show_spinner=False)
def load_local_sheet_preview_v2(path, mtime):
    df = pd.read_excel(path)
    return compute_match_preview(df)


def preview_button_label(title, preview=None, missing_text=None):
    if missing_text:
        return f"{title}\n:orange[{missing_text}]"
    if not preview:
        return f"{title}\n:orange[Préanalyse indisponible]"

    parts = []
    for system, label in (("iscn", "ISCN"), ("jon", "Jon")):
        stats = preview.get(system, {})
        if stats.get("percent") is not None:
            color = "green" if stats["percent"] == 100 else "red"
            parts.append(
                f":{color}[{label} {stats['ok']}/{stats['total']} "
                f"({format_percent(stats['percent'])}%)]"
            )
    if preview.get("errors"):
        parts.append(f":red[{preview['errors']} erreur(s)]")
    if not parts:
        parts.append(":orange[Sans références]")

    return f"{title}\n{' · '.join(parts)}"


def normalize_rule_catalog_frame(df):
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
            column: aliases.get(str(column).strip().lower(), str(column).strip())
            for column in df.columns
        }
    )


def normalized_catalog_text(value):
    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return ""
    return str(value).strip()


def compare_rule_catalogs(reference_df, distant_df):
    """Compare les lignes et les textes éditables des deux référentiels."""

    reference_df = normalize_rule_catalog_frame(reference_df)
    distant_df = normalize_rule_catalog_frame(distant_df)
    if "Rule ID" not in distant_df.columns:
        raise ValueError("Le fichier distant ne contient pas la colonne 'Rule ID'.")

    reference_rows = {
        normalized_catalog_text(row["Rule ID"]): row
        for _, row in reference_df.iterrows()
        if normalized_catalog_text(row.get("Rule ID"))
    }
    distant_rows = {
        normalized_catalog_text(row["Rule ID"]): row
        for _, row in distant_df.iterrows()
        if normalized_catalog_text(row.get("Rule ID"))
    }

    reference_ids = set(reference_rows)
    distant_ids = set(distant_rows)
    missing = sorted(reference_ids - distant_ids)
    extra = sorted(distant_ids - reference_ids)
    differences = []
    for rule_id in sorted(reference_ids & distant_ids):
        for column in ("Libellé", "Explication"):
            internal_value = normalized_catalog_text(reference_rows[rule_id].get(column))
            distant_value = normalized_catalog_text(distant_rows[rule_id].get(column))
            if internal_value != distant_value:
                differences.append(
                    {
                        "rule_id": rule_id,
                        "column": column,
                        "internal": internal_value,
                        "distant": distant_value,
                    }
                )
    return missing, extra, differences


def render_rule_catalog_table(reference_df, distant_df=None):
    """Affiche le référentiel et ses discordances dans un tableau unique."""

    columns = [
        "Rule ID",
        "Référentiel",
        "Score par défaut",
        "Critère technique",
        "Libellé",
        "Explication",
    ]
    reference_df = normalize_rule_catalog_frame(reference_df)
    has_distant_catalog = distant_df is not None
    distant_df = (
        normalize_rule_catalog_frame(distant_df)
        if has_distant_catalog
        else pd.DataFrame(columns=columns)
    )
    if has_distant_catalog:
        missing, extra, differences = compare_rule_catalogs(reference_df, distant_df)
    else:
        missing, extra, differences = [], [], []
    missing_set = set(missing)
    extra_set = set(extra)
    difference_keys = {
        (difference["rule_id"], difference["column"]): difference
        for difference in differences
    }

    reference_rows = {
        normalized_catalog_text(row.get("Rule ID")): row.to_dict()
        for _, row in reference_df.iterrows()
    }
    distant_rows = {
        normalized_catalog_text(row.get("Rule ID")): row.to_dict()
        for _, row in distant_df.iterrows()
        if normalized_catalog_text(row.get("Rule ID"))
    }
    reference_order = list(reference_rows)
    system_counters = {"ISCN": 0, "JON": 0}
    order_numbers = {}
    for rule_id in reference_order:
        system = normalized_catalog_text(reference_rows[rule_id].get("Référentiel"))
        suffix = "JON" if system == "Jondreville 2020" else "ISCN"
        system_counters[suffix] += 1
        order_numbers[rule_id] = (
            normalized_catalog_text(reference_rows[rule_id].get("N°"))
            or f"{system_counters[suffix]}_{suffix}"
        )

    # Les règles distantes inconnues sont insérées à leur position relative dans
    # le Google Sheet, avant la prochaine règle canonique connue.
    merged_ids = reference_order.copy()
    previous_known = None
    for rule_id in distant_rows:
        if rule_id in reference_rows:
            previous_known = rule_id
            continue
        if rule_id not in extra_set:
            continue
        insert_at = merged_ids.index(previous_known) + 1 if previous_known in merged_ids else 0
        while insert_at < len(merged_ids) and merged_ids[insert_at] in extra_set:
            insert_at += 1
        merged_ids.insert(insert_at, rule_id)
        previous_known = rule_id

    body_rows = []
    for rule_id in merged_ids:
        is_missing = rule_id in missing_set
        is_extra = rule_id in extra_set
        row_class = "catalog-row-missing" if is_missing else "catalog-row-extra" if is_extra else ""
        source = distant_rows.get(rule_id) if is_extra else reference_rows.get(rule_id, {})
        displayed = dict(source or {})
        if not is_missing and not is_extra:
            distant_row = distant_rows.get(rule_id, {})
            for column in ("Libellé", "Explication"):
                if column in distant_row:
                    displayed[column] = normalized_catalog_text(distant_row.get(column))

        displayed_system = normalized_catalog_text(displayed.get("Référentiel"))
        order_suffix = "JON" if displayed_system == "Jondreville 2020" else "ISCN"
        order_label = order_numbers.get(rule_id, f"—_{order_suffix}")
        cells = [
            f"<td class='catalog-rule-id'>{html.escape(rule_id)}</td>",
            f"<td class='catalog-order'>{order_label}</td>",
        ]
        for column in columns[1:]:
            value = normalized_catalog_text(displayed.get(column)) or "—"
            difference = difference_keys.get((rule_id, column))
            if difference:
                cells.append(
                    "<td class='catalog-diff-cell'>"
                    "<div class='catalog-version catalog-version-distant'>"
                    f"<strong>Distant :</strong> {html.escape(value)}</div>"
                    "<div class='catalog-version catalog-version-internal'>"
                    f"<strong>Interne :</strong> {html.escape(difference['internal']) or '—'}</div>"
                    "<div class='catalog-harmonize'>Harmoniser les deux fichiers.</div>"
                    "</td>"
                )
            else:
                cells.append(f"<td>{html.escape(value)}</td>")

        if is_missing:
            status = "Absente du Google Sheet : règle à ajouter dans le fichier distant."
            cells[-1] = cells[-1][:-5] + (
                "<div class='catalog-row-note'>" + html.escape(status) + "</div></td>"
            )
        elif is_extra:
            status = "Absente du CSV interne : règle à supprimer du fichier distant."
            cells[-1] = cells[-1][:-5] + (
                "<div class='catalog-row-note'>" + html.escape(status) + "</div></td>"
            )
        body_rows.append(f"<tr class='{row_class}'>{''.join(cells)}</tr>")

    headers = "".join(
        f"<th>{html.escape(column)}</th>"
        for column in ["Rule ID", "N°"] + columns[1:]
    )
    st.markdown(
        "<div class='catalog-table-scroll'><table class='catalog-main-table'>"
        f"<thead><tr>{headers}</tr></thead><tbody>{''.join(body_rows)}</tbody>"
        "</table></div>",
        unsafe_allow_html=True,
    )


# Titre de l'application
st.title("Analyseur de Formules Caryotypiques (ISCN)")

# Fonction pour créer un lien de téléchargement Excel
def get_excel_download_link(df, filename="resultats_analyse.xlsx"):
    """
    Crée un lien HTML pour télécharger un DataFrame en Excel
    """
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        df.to_excel(writer, index=False, sheet_name='Résultats')
    excel_data = output.getvalue()
    b64 = base64.b64encode(excel_data).decode()
    href = f'<a href="data:application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;base64,{b64}" download="{filename}" class="download-button">Télécharger les résultats en Excel</a>'
    return href

# Fonction pour formater les explications avec des puces colorées
def format_anomalies_html(anomalies_df):
    """
    Formate les anomalies avec des puces colorées pour l'affichage HTML
    """
    html = ""
    for _, row in anomalies_df.iterrows():
        anom = html.escape(str(row['Anomalie']))
        type_anom = html.escape(str(row['Type']))
        score = row['Score ISCN 2024']
        clones_raw = row['Clones']
        clones = html.escape(clones_raw.strip()) if isinstance(clones_raw, str) and clones_raw.strip() else '—'
        explication_raw = row['Explication']
        if isinstance(explication_raw, str) and explication_raw.strip():
            explication = html.escape(explication_raw.strip())
        elif explication_raw is None or (isinstance(explication_raw, float) and pd.isna(explication_raw)):
            explication = '—'
        else:
            explication = html.escape(str(explication_raw))
        
        # Couleur de la puce selon le type d'anomalie
        if score == 2:
            color = "#FF5733"  # Rouge pour les anomalies à 2 points
        elif score == 1:
            color = "#33A1FF"  # Bleu pour les anomalies à 1 point
        else:
            color = "#AAAAAA"  # Gris pour les anomalies à 0 point
        
        # Couleur du score
        score_color = "#FFFFFF"
        score_bg = "#555555"
        
        html += f"""
        <div style="margin-bottom: 10px; padding: 8px; border-left: 4px solid {color}; background-color: #f9f9f9;">
            <div style="display: flex; align-items: center; margin-bottom: 5px;">
                <span style="font-weight: bold; flex: 1;">{anom}</span>
                <span style="background-color: {score_bg}; color: {score_color}; border-radius: 12px; padding: 2px 8px; 
                      display: inline-block; font-weight: bold;">{score} pts</span>
            </div>
            <div style="margin-left: 10px; color: #666;">
                <div><strong>Type:</strong> {type_anom}</div>
                <div><strong>Clones:</strong> {clones}</div>
                <div><strong>Explication:</strong> {explication}</div>
            </div>
        </div>
        """
    return html

# Fonction pour un affichage compact similaire à l'analyse par fichier
def format_anomalies_compact(anomalies_df):
    """Affiche les anomalies regroupées par clone avec détails."""

    clone_blocks = {}
    clone_order = []

    def add_line(clone_label, line_html):
        key = clone_label or ""
        if key not in clone_blocks:
            clone_blocks[key] = []
            clone_order.append(key)
        clone_blocks[key].append(line_html)

    def clean_text(value):
        if isinstance(value, str) and value.strip():
            return html.escape(value.strip())
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return '—'
        return html.escape(str(value))

    def line_color(score):
        try:
            value = float(score)
        except (TypeError, ValueError):
            value = 0
        if value >= 1.5:
            return "#D75C37"
        if value >= 0.5:
            return "#3A6EA5"
        return "#9EA7B8"

    def build_rule_help(rule_id, applied_explanation):
        path = get_rule_path(str(rule_id or ""), public_text_url=RULE_CATALOG_SHEET_URL or None)
        if not path:
            return ""

        applied_text = ""
        if isinstance(applied_explanation, str) and applied_explanation.strip():
            applied_text = applied_explanation.strip()
        elif applied_explanation is not None and not (
            isinstance(applied_explanation, float) and pd.isna(applied_explanation)
        ):
            applied_text = str(applied_explanation)

        popover_id = f"rule-popover-{uuid.uuid4().hex}"
        items = []
        for step in path:
            selected = bool(step["selected"])
            step_class = " selected" if selected else ""
            marker = "retenue" if selected else "testée avant"
            score_text = step.get("default_score")
            score_html = f" · score {html.escape(str(score_text))}" if score_text != "" else ""
            open_attr = " open" if selected else ""
            technical_check = str(step.get("technical_check") or "").strip()
            technical_html = (
                f'<div class="rule-detail-block"><strong>Critère technique :</strong> {html.escape(technical_check)}</div>'
                if technical_check
                else ""
            )
            explanation = str(step.get("explanation") or "").strip()
            explanation_html = (
                f'<div class="rule-detail-block"><strong>Explication clinique :</strong> {html.escape(explanation)}</div>'
                if explanation
                else ""
            )
            items.append(
                f'<li class="rule-step{step_class}">'
                f'<span class="rule-step-order">{step["order"]}</span>'
                f'<details class="rule-details"{open_attr}>'
                f'<summary>'
                f'<span><strong>{html.escape(str(step["rule_id"]))}</strong>'
                f' <em>{html.escape(marker)}</em>{score_html}<br>'
                f'<span class="rule-title">{html.escape(str(step["title"]))}</span></span>'
                f'</summary>'
                f'<div class="rule-detail-panel">'
                f'{explanation_html}'
                f'{technical_html}'
                f'</div>'
                f'</details>'
                f'</li>'
            )

        applied_html = ""
        if applied_text:
            applied_html = (
                '<div class="rule-applied">'
                f'<strong>Explication appliquée :</strong> {html.escape(applied_text)}'
                '</div>'
            )

        return (
            f'<button type="button" class="rule-help" popovertarget="{popover_id}" '
            'aria-label="Afficher le parcours de règles">?</button>'
            f'<div id="{popover_id}" class="rule-popover" popover>'
            f'{applied_html}'
            '<ol>'
            f'{"".join(items)}'
            '</ol>'
            '</div>'
        )

    def build_pill_text(label, score, explanation, detail_text, rule_id, applied_explanation):
        detail = detail_text
        if explanation != '—':
            detail = f"{detail_text} ({explanation})"
        rule_help = build_rule_help(rule_id, applied_explanation)
        score_value = render_score_value(score)
        return f"{label} {rule_help}{score_value} : {detail}"

    for _, row in anomalies_df.iterrows():
        score_iscn = row['Score ISCN 2024']
        score_jon = row['Score Jondreville 2020']
        anomalie = html.escape(str(row['Anomalie']))
        type_text = clean_text(row.get('Type'))
        if type_text == '—':
            type_text = anomalie

        explication_iscn = clean_text(row.get('Explication'))
        explication_jon = clean_text(row.get('Explication Jondreville 2020'))
        rule_id_iscn = row.get('RuleID_ISCN') or ''
        rule_id_jon = row.get('RuleID_Jon') or ''
        rule_explanation_iscn = row.get('RuleExplanation_ISCN') or row.get('Explication')
        rule_explanation_jon = row.get('RuleExplanation_Jon') or row.get('Explication Jondreville 2020')
        clone_details = row.get('CloneDetails')
        if not isinstance(clone_details, list) or not clone_details:
            clone_details = [{"label": None, "is_reference": True, "reason": "", "count": 1}]

        for detail in clone_details:
            label = detail.get('label')
            is_reference = detail.get('is_reference', True)
            reason_raw = detail.get('reason') or ''
            reason_html = clean_text(reason_raw) if reason_raw else ''
            count = detail.get('count', 1)

            score_override = detail.get('score_override')
            if score_override is not None:
                line_score_iscn = score_override
                line_score_jon = score_override
                line_exp_iscn = reason_html or explication_iscn
                line_exp_jon = reason_html or explication_jon
            elif is_reference:
                line_score_iscn = score_iscn
                line_exp_iscn = explication_iscn
                line_score_jon = score_jon
                line_exp_jon = explication_jon
            else:
                line_score_iscn = 0
                line_exp_iscn = reason_html or explication_iscn
                line_score_jon = 0
                line_exp_jon = reason_html or explication_jon

            color = line_color(line_score_iscn)
            anomaly_label = f"[{anomalie}]"
            if count > 1:
                anomaly_label += f" &times;{count}"
            line_html = (
                f'<div class="anomaly-line" style="border-left-color: {color};">'
                f'<span class="anomaly-label" style="color: {color};">{anomaly_label}</span>'
                f'<span class="score-pill score-pill-iscn">{build_pill_text("ISCN", line_score_iscn, line_exp_iscn, type_text, rule_id_iscn, rule_explanation_iscn)}</span>'
                f'<span class="score-pill score-pill-jon">{build_pill_text("Jon", line_score_jon, line_exp_jon, type_text, rule_id_jon, rule_explanation_jon)}</span>'
                '</div>'
            )

            add_line(label, line_html)

    if not clone_blocks:
        return ""

    labeled_keys = [lbl for lbl in clone_order if lbl]
    show_labels = len(labeled_keys) >= 2

    blocks = []
    for label in clone_order:
        lines_html = "".join(clone_blocks[label])
        if show_labels and label:
            blocks.append(
                '<div class="anomaly-compact">'
                f'<span class="clone-pill">{html.escape(label)}</span>'
                f'<div class="clone-group-lines">{lines_html}</div>'
                '</div>'
            )
        else:
            blocks.append(f'<div class="anomaly-compact">{lines_html}</div>')

    return "\n".join(blocks)


def score_tone_class(score):
    try:
        value = float(score)
    except (TypeError, ValueError):
        return "score-tone-default"
    if value >= 1.5:
        return "score-tone-2"
    if value >= 0.5:
        return "score-tone-1"
    return "score-tone-0"


def format_score_text(score):
    if isinstance(score, str):
        return score
    if score is None or (isinstance(score, float) and pd.isna(score)):
        return '—'
    if isinstance(score, float) and score.is_integer():
        return str(int(score))
    return str(score)


def render_score_value(score):
    text = html.escape(format_score_text(score))
    tone_class = score_tone_class(score)
    return f"<span class='score-value {tone_class}'>{text}</span>"


def render_plain_score(score):
    text = html.escape(format_score_text(score))
    return f"<span class='score-value score-tone-default'>{text}</span>"


def render_score_totals(score_iscn, score_jon):
    """Crée le bloc HTML présentant les scores globaux."""
    return f"""
    <div class="score-summary-group">
        <div class="score-summary">
            <span class="score-label">ISCN</span>
            <span class="score-pill score-pill-iscn">{render_plain_score(score_iscn)}</span>
        </div>
        <div class="score-summary">
            <span class="score-label">Jondreville</span>
            <span class="score-pill score-pill-jon">{render_plain_score(score_jon)}</span>
        </div>
    </div>
    """

# Interface utilisateur
st.markdown("""
Cette application permet d'analyser des formules caryotypiques (notation ISCN) pour :
- Compter le nombre d'anomalies
- Identifier le type de chaque anomalie
- Comparer le comptage automatique avec un comptage manuel (si disponible)
""")

with st.expander("Référentiel des règles de scoring"):
    if RULE_CATALOG_SHEET_URL:
        st.markdown(
            f'<a href="{html.escape(RULE_CATALOG_SHEET_URL)}" target="_blank" '
            'rel="noopener noreferrer">Ouvrir le Google Sheet du catalogue</a>',
            unsafe_allow_html=True,
        )
    try:
        internal_catalog = pd.read_csv(RULE_CATALOG_REFERENCE_PATH)
    except Exception as catalog_error:
        internal_catalog = get_canonical_rule_catalog_dataframe()
        st.error(f"CSV interne du référentiel indisponible : {catalog_error}")

    # Le CSV doit évoluer dans le même commit que RULE_CATALOG.
    code_catalog = get_canonical_rule_catalog_dataframe()
    snapshot_columns = [
        "Rule ID",
        "N°",
        "Référentiel",
        "Score par défaut",
        "Critère technique",
        "Libellé",
        "Explication",
    ]
    internal_snapshot = internal_catalog[snapshot_columns].fillna("").astype(str)
    code_snapshot = code_catalog[snapshot_columns].fillna("").astype(str)
    if internal_snapshot.to_dict("records") != code_snapshot.to_dict("records"):
        st.error(
            "Le code et rules_catalog_reference.csv ne sont plus synchronisés. "
            "Exécuter : python scripts/sync_rule_catalog_reference.py"
        )

    distant_catalog, distant_catalog_error = load_google_sheet(RULE_CATALOG_SHEET_URL)
    if distant_catalog_error:
        st.warning(
            "Comparaison avec le fichier distant impossible : "
            f"{distant_catalog_error}"
        )
        render_rule_catalog_table(internal_catalog)
    else:
        try:
            render_rule_catalog_table(internal_catalog, distant_catalog)
        except Exception as comparison_error:
            st.warning(f"Contrôle du référentiel distant impossible : {comparison_error}")
            render_rule_catalog_table(internal_catalog)

# Création des onglets (par défaut: analyse d'un fichier)
tab2, tab1 = st.tabs(["Analyse d'un fichier", "Analyse d'une formule"])

# Onglet 1: Analyse d'une formule
with tab1:
    st.subheader("Entrez une formule caryotypique")
    formule = st.text_input("Formule ISCN", placeholder="Ex: 47,XX,+8[20]")
    
    if st.button("Analyser la formule", key="analyser_formule"):
        if formule:
            df, totals, error = analyser_formule(formule, debug=True)
            if error:
                st.error(error)
            else:
                st.success(
                    f"Scores détectés — ISCN: {totals['iscn']} | Jondreville: {totals['jondroville']}"
                )

                st.markdown(
                    render_score_totals(totals['iscn'], totals['jondroville']),
                    unsafe_allow_html=True
                )

                if totals.get("formule_equivalente"):
                    st.markdown(
                        f"""
                        <div style="margin: 10px 0;">
                            <div><strong>Formule :</strong> {html.escape(totals.get('formule_originale', ''))}<br>
                                <span style="display:inline-block; border:2px solid #f5c542; border-radius:999px; padding:2px 8px; background:#fff8e1;">
                                    {html.escape(totals.get('formule_equivalente', ''))}
                                </span>
                            </div>
                        </div>
                        """,
                        unsafe_allow_html=True,
                    )

                # Affichage du tableau avec info-bulles
                st.markdown("### Détail des anomalies")

                # Formatage des anomalies pour l'affichage
                anomalies_df = df.iloc[:-1]  # Exclure la ligne TOTAL
                anomalies_html = format_anomalies_compact(anomalies_df)
                st.markdown(anomalies_html, unsafe_allow_html=True)

        else:
            st.warning("Veuillez entrer une formule caryotypique.")

# Onglet 2: Analyse d'un fichier
with tab2:
    st.subheader("Chargez un fichier contenant des formules caryotypiques")
    st.markdown(
        f'<a href="{html.escape(TEST_SHEET_URL)}" target="_blank" '
        'rel="noopener noreferrer">Ouvrir le Google Sheet des formules de test</a>',
        unsafe_allow_html=True,
    )
    remote_preview, remote_preview_err = load_google_sheet_preview_v2(TEST_SHEET_URL)
    remote_button_label = preview_button_label(
        "Analyser le fichier de tests",
        remote_preview,
        "Préanalyse indisponible" if remote_preview_err else None,
    )
    if LOCAL_TEST_PATH.exists():
        local_preview = load_local_sheet_preview_v2(
            str(LOCAL_TEST_PATH), LOCAL_TEST_PATH.stat().st_mtime
        )
        local_button_label = preview_button_label("Analyser le fichier local MYC", local_preview)
    else:
        local_button_label = preview_button_label(
            "Analyser le fichier local MYC",
            missing_text="En attente du fichier MYC",
        )

    test_col, local_col = st.columns([1, 1], gap="small")
    test_button = test_col.button(remote_button_label, width="stretch")
    local_test_button = local_col.button(local_button_label, width="stretch")

    show_local_uploader = st.session_state.get("show_local_myc_uploader", False)
    if show_local_uploader:
        st.warning(
            "Pour des raisons de confidentialité, le fichier local MYC n'est pas stocké en ligne. "
            f"Placez `{LOCAL_TEST_FILENAME}` dans le répertoire de l'application ou chargez-le ici. "
            "Le fichier chargé est lu en mémoire pour cette session et n'est pas écrit sur disque par l'application."
        )
        local_uploaded_file = st.file_uploader(
            "Charger le fichier local MYC", type=["csv", "xlsx", "xls"], key="local_myc_file"
        )
    else:
        local_uploaded_file = None

    uploaded_file = st.file_uploader(
        "Choisir un autre fichier CSV ou Excel", type=["csv", "xlsx", "xls"], key="file"
    )
 

    # Préserver les données chargées entre les reruns (ex: clic sur "Trier")
    df_input = st.session_state.get("df_input")
    if test_button:
        st.session_state["show_local_myc_uploader"] = False
        df_input, err = load_google_sheet(TEST_SHEET_URL)
        if err:
            st.error(f"Erreur lors du chargement du fichier de tests : {err}")
            st.stop()
        else:
            st.session_state["df_input"] = df_input
    elif local_test_button:
        if LOCAL_TEST_PATH.exists():
            df_input = pd.read_excel(LOCAL_TEST_PATH)
            st.session_state["df_input"] = df_input
            st.session_state["show_local_myc_uploader"] = False
        else:
            st.session_state["show_local_myc_uploader"] = True
            st.session_state["df_input"] = None
            st.warning(
                "Pour des raisons de confidentialité, le fichier n'est pas stocké en ligne. "
                f"Veuillez charger votre fichier de test ou déposer `{LOCAL_TEST_FILENAME}` dans le répertoire courant."
            )
            st.rerun()
    elif local_uploaded_file is not None:
        if local_uploaded_file.name.endswith('.csv'):
            df_input = pd.read_csv(local_uploaded_file)
        else:
            df_input = pd.read_excel(local_uploaded_file)
        st.session_state["df_input"] = df_input
        st.session_state["show_local_myc_uploader"] = False
    elif uploaded_file is not None:
        st.session_state["show_local_myc_uploader"] = False
        # Déterminer le type de fichier
        if uploaded_file.name.endswith('.csv'):
            df_input = pd.read_csv(uploaded_file)
        else:  # Excel
            df_input = pd.read_excel(uploaded_file)
        st.session_state["df_input"] = df_input



    if df_input is not None:
        try:
            
            # Rechercher la colonne Formule de manière insensible à la casse
            formule_col = None
            for col in df_input.columns:
                if col.strip().lower() == 'formule':
                    formule_col = col
                    break

            if formule_col is None:
                st.error("Le fichier doit contenir au moins une colonne 'Formule'.")
            else:
                # Renommer la colonne trouvée en 'Formule' pour simplifier la suite
                if formule_col != 'Formule':
                    df_input = df_input.rename(columns={formule_col: 'Formule'})

                # Détection et renommage des colonnes de référence
                count_i_col = None
                count_j_col = None
                for col in df_input.columns:
                    norm = col.strip().lower()
                    if norm in ('count_i', 'count_iscn'):
                        count_i_col = col
                    elif norm in ('count_j', 'count_jon'):
                        count_j_col = col

                if count_i_col and count_i_col != 'Count_i':
                    df_input = df_input.rename(columns={count_i_col: 'Count_i'})
                if count_j_col and count_j_col != 'Count_j':
                    df_input = df_input.rename(columns={count_j_col: 'Count_j'})

                has_count_i = 'Count_i' in df_input.columns
                has_count_j = 'Count_j' in df_input.columns

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

                # Création du DataFrame de résultats
                results = []
                all_anomalies_details = []
                match_details = []

                for idx, row in df_input.iterrows():
                    formule_fichier = row['Formule']
                    ref_iscn_value = normalize_reference(row['Count_i']) if has_count_i else None
                    ref_jondroville_value = normalize_reference(row['Count_j']) if has_count_j else None

                    df_analyse, totals, error = analyser_formule(formule_fichier, debug=True)

                    match_detail = {"iscn": None, "jon": None}

                    if error:
                        anomalies_detail = error
                        comptage_iscn = "Erreur"
                        comptage_jondroville = "Erreur"
                        all_anomalies_details.append({"error": True, "message": error})
                    else:
                        # Extraction des détails des anomalies
                        anomalies_df = df_analyse.iloc[:-1]  # Exclure la ligne TOTAL

                        # Stocker les détails pour l'affichage
                        all_anomalies_details.append({"error": False, "df": anomalies_df})

                        # Texte simple pour l'export
                        def clean_detail(value):
                            if isinstance(value, str) and value.strip():
                                return value.strip()
                            if value is None or (isinstance(value, float) and pd.isna(value)):
                                return "—"
                            return str(value)

                        detail_chunks = []
                        for _, row_detail in anomalies_df.iterrows():
                            type_label = clean_detail(row_detail.get('Type'))
                            if type_label == '—':
                                type_label = clean_detail(row_detail.get('Anomalie'))
                            exp_iscn = clean_detail(row_detail.get('Explication'))
                            exp_jon = clean_detail(row_detail.get('Explication Jondreville 2020'))
                            clone_details = row_detail.get('CloneDetails')
                            if not isinstance(clone_details, list) or not clone_details:
                                clone_details = [{"label": None, "is_reference": True, "reason": "", "count": 1}]

                            multiple_clones = len([cd for cd in clone_details if cd.get('label')]) >= 2

                            for detail in clone_details:
                                is_reference = detail.get('is_reference', True)
                                reason = clean_detail(detail.get('reason')) if detail.get('reason') else ''
                                clone_label = detail.get('label') if multiple_clones else ''
                                count = detail.get('count', 1)

                                score_override = detail.get('score_override')
                                if score_override is not None:
                                    line_score_iscn = score_override
                                    line_score_jon = score_override
                                    line_exp_iscn = reason or exp_iscn
                                    line_exp_jon = reason or exp_jon
                                else:
                                    line_score_iscn = row_detail['Score ISCN 2024'] if is_reference else 0
                                    line_exp_iscn = exp_iscn if is_reference else (reason or exp_iscn)
                                    line_score_jon = row_detail['Score Jondreville 2020'] if is_reference else 0
                                    line_exp_jon = exp_jon if is_reference else (reason or exp_jon)

                                prefix = f"{clone_label}: " if clone_label else ""
                                label_with_count = f"{type_label} x{count}" if count > 1 else type_label
                                chunk = (
                                    f"{prefix}{label_with_count}: ISCN {line_score_iscn}"
                                    f" ({line_exp_iscn}) | Jon {line_score_jon}"
                                    f" ({line_exp_jon})"
                                )
                                detail_chunks.append(chunk)

                        anomalies_detail = ", ".join(detail_chunks)

                        comptage_iscn = totals['iscn']
                        comptage_jondroville = totals['jondroville']

                        if has_count_i and ref_iscn_value is not None:
                            match_detail['iscn'] = comptage_iscn == ref_iscn_value
                        if has_count_j and ref_jondroville_value is not None:
                            match_detail['jon'] = comptage_jondroville == ref_jondroville_value

                    match_details.append(match_detail)

                    result_row = {
                        "Ligne": idx + 1,  # +1 pour ne pas inclure l'en-tête du fichier
                        "Formule": formule_fichier,
                        "Comptage ISCN": comptage_iscn,
                        "Comptage Jon": comptage_jondroville,
                        "Anomalies détectées": anomalies_detail  # Version texte pour l'export
                    }
                    if totals.get("formule_equivalente"):
                        result_row["_FormuleEquivalente"] = totals.get("formule_equivalente", "")

                    if has_count_i:
                        result_row["Ref ISCN"] = ref_iscn_value
                    if has_count_j:
                        result_row["Ref Jon"] = ref_jondroville_value

                    results.append(result_row)

                # Création du DataFrame de résultats avec un ordre de colonnes défini
                results_df = pd.DataFrame(results)
                columns_order = ["Ligne", "Formule", "Comptage ISCN"]
                if has_count_i:
                    columns_order.append("Ref ISCN")
                columns_order.append("Comptage Jon")
                if has_count_j:
                    columns_order.append("Ref Jon")
                if "_FormuleEquivalente" in results_df.columns:
                    columns_order.append("_FormuleEquivalente")
                columns_order.append("Anomalies détectées")
                results_df = results_df[columns_order]

                # Fonctions utilitaires pour l'affichage
                def format_display(value):
                    if isinstance(value, str):
                        return value
                    if pd.isna(value):
                        return "—"
                    if isinstance(value, float) and value.is_integer():
                        return str(int(value))
                    return str(value)

                display_labels = ["Ligne", "Formule", "Comptage ISCN"]
                if has_count_i:
                    display_labels.append("Ref ISCN")
                display_labels.append("Comptage Jon")
                if has_count_j:
                    display_labels.append("Ref Jon")
                display_labels.append("Anomalies détectées")

                # Affichage des résultats
                title_col, sort_col = st.columns([1, 0.2])
                title_col.markdown("### Résultats de l'analyse")

                # Option de tri par discordance, uniquement si au moins une référence est présente
                sort_enabled = False
                if has_count_i or has_count_j:
                    sort_enabled = sort_col.checkbox(
                        "Trier",
                        value=False,
                        help="Affiche d'abord les lignes discordantes (ISCN ou Jon), dans l'ordre des lignes d'origine."
                    )

                # Longueurs de sécurité pour éviter toute désynchronisation
                base_len = min(len(results_df), len(match_details), len(all_anomalies_details))
                display_order = list(range(base_len))
                if sort_enabled:
                    if base_len < len(results_df) or base_len < len(match_details):
                        st.warning("Tri non appliqué: longueurs incohérentes des données d'affichage.")
                    else:
                        try:
                            def sort_key(idx: int):
                                ligne_val = results_df.at[idx, "Ligne"] if "Ligne" in results_df else idx
                                is_discordant = (
                                    match_details[idx].get("iscn") is False
                                    or match_details[idx].get("jon") is False
                                )
                                return (0 if is_discordant else 1, ligne_val)

                            display_order = sorted(display_order, key=sort_key)
                        except Exception as sort_err:
                            st.warning(f"Tri non appliqué (erreur: {sort_err})")

                header_cells = []
                for label in display_labels:
                    if label == "Ligne":
                        header_cells.append(
                            '<th class="line-jump-header">'
                            '<div class="line-jump-title">Ligne</div>'
                            '<input class="line-jump-input" type="number" min="1" '
                            'placeholder="N°" aria-label="Aller à une ligne">'
                            '</th>'
                        )
                    else:
                        header_cells.append(f"<th>{html.escape(label)}</th>")
                header_html = "".join(header_cells)

                body_rows = []
                for idx in display_order:
                    row_data = results_df.iloc[idx]
                    line_number = html.escape(format_display(row_data.get("Ligne")))
                    anomalies = all_anomalies_details[idx]
                    matches = match_details[idx] if idx < len(match_details) else {"iscn": None, "jon": None}
                    cells = []

                    for label in display_labels:
                        if label == "Ligne":
                            cells.append(f'<td class="line-number-cell">{line_number}</td>')
                        elif label == "Comptage ISCN":
                            raw_value = row_data.get(label)
                            icon = ""
                            if matches["iscn"] is not None:
                                icon = f"<span class='pill-icon'>{'✅' if matches['iscn'] else '❌'}</span>"
                            badge = f"<span class='score-pill score-pill-iscn'>{render_plain_score(raw_value)}{icon}</span>"
                            cells.append(f"<td class='score-cell'>{badge}</td>")
                        elif label == "Comptage Jon":
                            raw_value = row_data.get(label)
                            icon = ""
                            if matches["jon"] is not None:
                                icon = f"<span class='pill-icon'>{'✅' if matches['jon'] else '❌'}</span>"
                            badge = f"<span class='score-pill score-pill-jon'>{render_plain_score(raw_value)}{icon}</span>"
                            cells.append(f"<td class='score-cell'>{badge}</td>")
                        elif label == "Anomalies détectées":
                            if anomalies["error"]:
                                message = html.escape(anomalies["message"])
                                cells.append(
                                    f"<td><div class=\"anomaly-error\">{message}</div></td>"
                                )
                            else:
                                anomalies_html = format_anomalies_compact(anomalies["df"])
                                cells.append(f"<td>{anomalies_html}</td>")
                        elif label == "Formule":
                            value = format_display(row_data.get(label))
                            equiv = format_display(row_data.get("_FormuleEquivalente"))
                            if equiv and equiv != "—":
                                ring = (
                                    "<span style='display:inline-block; border:2px solid #f5c542; "
                                    "border-radius:999px; padding:2px 8px; background:#fff8e1;'>"
                                    f"{html.escape(equiv)}</span>"
                                )
                                cells.append(f"<td>{html.escape(value)}<br>{ring}</td>")
                            else:
                                cells.append(f"<td>{html.escape(value)}</td>")
                        else:
                            value = format_display(row_data.get(label))
                            cells.append(f"<td>{html.escape(value)}</td>")

                    body_rows.append(f'<tr id="formule-{line_number}">{"".join(cells)}</tr>')

                table_html = f"""
                <div class="results-table-container">
                    <table class="results-table">
                        <thead>
                            <tr>{header_html}</tr>
                        </thead>
                        <tbody>
                            {''.join(body_rows)}
                        </tbody>
                    </table>
                </div>
                """

                st.markdown(table_html, unsafe_allow_html=True)
                components.html(
                    """
                    <script>
                    (function () {
                        function bindLineJumpInput() {
                            const doc = window.parent.document;
                            const inputs = doc.querySelectorAll("th.line-jump-header input.line-jump-input:not([data-line-jump-bound='1'])");
                            inputs.forEach(function (input) {
                                input.dataset.lineJumpBound = "1";
                                input.addEventListener("focus", function () {
                                    this.select();
                                });
                                input.addEventListener("keydown", function (event) {
                                    if (event.key !== "Enter") {
                                        return;
                                    }
                                    event.preventDefault();
                                    const value = this.value.trim();
                                    if (!/^\\d+$/.test(value)) {
                                        return;
                                    }
                                    const row = doc.getElementById("formule-" + value);
                                    if (!row) {
                                        return;
                                    }
                                    doc.querySelectorAll(".formula-row-target").forEach(function (activeRow) {
                                        activeRow.classList.remove("formula-row-target");
                                    });
                                    row.scrollIntoView({block: "center", behavior: "smooth"});
                                    row.classList.add("formula-row-target");
                                });
                            });
                        }
                        bindLineJumpInput();
                        window.setTimeout(bindLineJumpInput, 500);
                    })();
                    </script>
                    """,
                    height=0,
                )

                # Statistiques de correspondance + export sur une seule ligne
                columns = st.columns(3)
                renderers = []
                match_summary = summarize_matches(match_details)

                if has_count_i:
                    total_i = match_summary["iscn"]["total"]
                    match_i = match_summary["iscn"]["ok"]
                    if total_i:
                        percent_i = format_percent(match_summary["iscn"]["percent"])
                        msg_i = f"Correspondance ISCN: {match_i}/{total_i} ({percent_i}%)"

                        def render_iscn(col, message=msg_i):
                            col.success(message)

                        renderers.append(render_iscn)

                if has_count_j:
                    total_j = match_summary["jon"]["total"]
                    match_j = match_summary["jon"]["ok"]
                    if total_j:
                        percent_j = format_percent(match_summary["jon"]["percent"])
                        msg_j = f"Correspondance Jondreville: {match_j}/{total_j} ({percent_j}%)"

                        def render_jon(col, message=msg_j):
                            col.success(message)

                        renderers.append(render_jon)

                download_html = get_excel_download_link(results_df)
                download_html = download_html.replace(
                    'class="download-button">',
                    'class="download-button"><span class="icon">⬇️</span>'
                )

                def render_download(col, html=download_html):
                    col.markdown(html, unsafe_allow_html=True)

                renderers.append(render_download)

                for idx, render in enumerate(renderers[:3]):
                    render(columns[idx])

                for idx in range(len(renderers), 3):
                    columns[idx].markdown("&nbsp;", unsafe_allow_html=True)
                
        except Exception as e:
            st.error(f"Erreur lors de l'analyse du fichier: {str(e)}")

# CSS pour améliorer l'apparence
st.markdown("""
<style>
    .download-button {
        display: inline-flex;
        align-items: center;
        gap: 10px;
        padding: 12px 28px;
        background-color: #111827;
        color: #ffffff;
        text-decoration: none;
        border-radius: 999px;
        margin-top: 10px;
        font-weight: 600;
        text-align: center;
        transition: background-color 0.2s ease, transform 0.2s ease;
    }
    
    .download-button:hover {
        background-color: #000000;
        transform: translateY(-1px);
    }

    .download-button .icon {
        font-size: 1.1rem;
    }

    .score-summary-group {
        display: flex;
        gap: 16px;
        flex-wrap: wrap;
        margin: 10px 0 20px;
    }

    .score-summary {
        display: flex;
        align-items: center;
        gap: 8px;
        background-color: #f4f5f7;
        border-radius: 999px;
        padding: 6px 12px;
    }

    .score-label {
        font-weight: 600;
        color: #1f2937;
    }

    .score-pill {
        display: inline-flex;
        align-items: center;
        position: relative;
        border-radius: 999px;
        padding: 2px 10px;
        font-weight: 600;
        font-size: 0.9rem;
        gap: 4px;
        background-color: transparent;
        border: 1px solid currentColor;
    }

    .score-pill-iscn {
        border-color: #3a6ea5;
        color: #1f2937;
    }

    .score-pill-jon {
        border-color: #8a5a9e;
        color: #1f2937;
    }

    h3 {
        margin-top: 30px;
        margin-bottom: 20px;
        color: #1E3A8A;
    }

    .results-table-container {
        max-height: 900px;
        overflow-y: auto;
        overflow-x: auto;
        border: 1px solid #e5e7eb;
        border-radius: 6px;
        margin-top: 1rem;
        width: 100%;
        max-width: 100%;
    }

    .results-table {
        width: 100%;
        border-collapse: collapse;
        table-layout: auto;
    }

    .results-table th,
    .results-table td {
        padding: 12px;
        vertical-align: top;
        border-bottom: 1px solid #e5e7eb;
        background-color: #ffffff;
        overflow-wrap: anywhere;
        word-break: break-word;
    }

    .results-table tr:target td {
        background-color: #fff7d6;
        box-shadow: inset 4px 0 0 #f5c542;
    }

    .results-table tr.formula-row-target td {
        background-color: #fff7d6;
        box-shadow: inset 4px 0 0 #f5c542;
    }

    .line-number-cell {
        white-space: nowrap;
        font-weight: 700;
        color: #111827;
        text-align: center;
    }

    .line-jump-header {
        min-width: 86px;
    }

    .line-jump-title {
        margin-bottom: 4px;
    }

    .line-jump-input {
        width: 64px;
        height: 28px;
        border: 1px solid #cbd5e1;
        border-radius: 4px;
        padding: 2px 6px;
        font: inherit;
        font-weight: 700;
        text-align: center;
        background: #ffffff;
        color: #111827;
    }

    .line-jump-input:focus {
        border-color: #3a6ea5;
        box-shadow: 0 0 0 2px rgba(58, 110, 165, 0.16);
        outline: none;
    }

    /* Ajuster la largeur des colonnes Formule et Anomalies détectées */
    .results-table th:nth-child(2),
    .results-table td:nth-child(2) {
        width: 28%;
        max-width: 28%;
        white-space: normal;
    }

    .results-table th:last-child,
    .results-table td:last-child {
        width: 40%;
        min-width: 40%;
        white-space: normal;
    }

    .results-table th {
        position: sticky;
        top: 0;
        background-color: #f1f5ff;
        z-index: 2;
        text-align: left;
        font-weight: 600;
    }

    .results-table tbody tr:nth-child(even) td {
        background-color: #f9fafb;
    }

    .anomaly-error {
        color: #b91c1c;
        font-weight: 600;
    }

    .anomaly-compact {
        display: flex;
        flex-direction: column;
        align-items: flex-start;
        gap: 6px;
        margin: 4px 0;
        padding: 8px 10px;
        background-color: #f9f9f9;
        border-radius: 6px;
        font-size: 14px;
        line-height: 1.4;
    }

    .clone-group {
        display: flex;
        flex-direction: column;
        gap: 6px;
        width: 100%;
    }

    .clone-group-lines {
        display: flex;
        flex-direction: column;
        gap: 6px;
        width: 100%;
    }

    .anomaly-line {
        display: flex;
        flex-wrap: wrap;
        align-items: center;
        gap: 6px;
        width: 100%;
        background-color: #ffffff;
        border-radius: 4px;
        padding: 6px 8px;
        border-left: 3px solid transparent;
    }

    .clone-pill {
        background-color: #111827;
        color: #ffffff;
        border-radius: 999px;
        padding: 2px 12px;
        font-size: 12px;
        font-weight: 600;
        align-self: flex-start;
    }

    .anomaly-compact .anomaly-label {
        font-weight: 600;
    }

    .anomaly-compact .score-pill {
        font-size: 12px;
        padding: 2px 8px;
    }

    .pill-icon {
        font-size: 0.85em;
        margin-left: 4px;
    }

    .score-value {
        font-weight: 700;
    }

    .score-tone-2 {
        color: #D75C37;
    }

    .score-tone-1 {
        color: #3A6EA5;
    }

    .score-tone-0 {
        color: #9EA7B8;
    }

    .score-tone-default {
        color: #1f2937;
    }

    .score-cell {
        white-space: nowrap;
    }

    .score-cell .score-pill {
        font-size: 0.85rem;
    }

    .rule-help {
        width: 18px;
        height: 18px;
        display: inline-flex;
        align-items: center;
        justify-content: center;
        flex: 0 0 auto;
        margin: 0 2px;
        padding: 0;
        border-radius: 999px;
        background-color: #ffffff;
        color: #374151;
        border: 1px solid #9ca3af;
        font-size: 12px;
        font-weight: 800;
        cursor: pointer;
        line-height: 1;
        font-family: inherit;
    }

    .rule-help:hover {
        background-color: #f8fafc;
        border-color: #64748b;
    }

    .rule-popover {
        min-width: 360px;
        max-width: min(560px, 85vw);
        max-height: 420px;
        overflow: auto;
        padding: 12px;
        background-color: #ffffff;
        color: #111827;
        border: 1px solid #d1d5db;
        border-radius: 6px;
        box-shadow: 0 12px 30px rgba(15, 23, 42, 0.18);
        white-space: normal;
        font-weight: 400;
        margin: auto;
    }

    .rule-popover ol {
        margin: 8px 0 0;
        padding: 0;
        list-style: none;
    }

    .rule-step {
        display: grid;
        grid-template-columns: 24px 1fr;
        gap: 8px;
        padding: 6px 0;
        border-top: 1px solid #eef2f7;
    }

    .rule-step:first-child {
        border-top: none;
    }

    .rule-step.selected {
        color: #0f5132;
    }

    .rule-step:target {
        outline: 2px solid #2563eb;
        outline-offset: 2px;
        border-radius: 4px;
        background-color: #eff6ff;
    }

    .rule-step-order {
        width: 20px;
        height: 20px;
        display: inline-flex;
        align-items: center;
        justify-content: center;
        border-radius: 999px;
        background-color: #eef2f7;
        color: #475569;
        font-size: 11px;
        font-weight: 700;
    }

    .rule-step.selected .rule-step-order {
        background-color: #d1fae5;
        color: #065f46;
    }

    .rule-details {
        min-width: 0;
    }

    .rule-details summary {
        cursor: pointer;
        list-style-position: outside;
        color: #111827;
    }

    .rule-step.selected .rule-details summary {
        color: #0f5132;
        font-weight: 600;
    }

    .rule-title {
        color: #334155;
        font-size: 0.84rem;
    }

    .rule-detail-panel {
        margin: 6px 0 0 14px;
        padding: 7px 8px;
        border-left: 2px solid #cbd5e1;
        background-color: #f8fafc;
        border-radius: 4px;
    }

    .rule-step.selected .rule-detail-panel {
        border-left-color: #10b981;
        background-color: #ecfdf5;
    }

    .rule-detail-block {
        margin-top: 5px;
        color: #475569;
        font-size: 0.78rem;
        line-height: 1.38;
    }

    .rule-detail-block:first-child {
        margin-top: 0;
    }

    .rule-applied {
        padding: 8px;
        border-radius: 4px;
        background-color: #f8fafc;
        border: 1px solid #e2e8f0;
        font-size: 0.85rem;
    }

    .catalog-table-scroll {
        width: 100%;
        max-height: 70vh;
        overflow: auto;
        margin-top: 10px;
        border: 1px solid #d1d5db;
    }

    .catalog-main-table {
        width: 100%;
        margin-bottom: 12px;
        border-collapse: collapse;
        font-size: 0.82rem;
    }

    .catalog-main-table th,
    .catalog-main-table td {
        padding: 7px;
        border: 1px solid #d1d5db;
        vertical-align: top;
        text-align: left;
    }

    .catalog-main-table th {
        position: sticky;
        top: 0;
        z-index: 2;
        background: #f1f5f9;
        color: #0f172a;
    }

    .catalog-main-table th:first-child {
        left: 0;
        z-index: 4;
    }

    .catalog-rule-id {
        position: sticky;
        left: 0;
        z-index: 1;
        min-width: 245px;
        background: #ffffff;
        font-weight: 700;
    }

    .catalog-order {
        min-width: 38px;
        text-align: center !important;
        font-weight: 800;
    }

    .catalog-row-missing td {
        background-color: #dbeafe;
        color: #1e3a8a;
    }

    .catalog-row-missing .catalog-rule-id {
        background-color: #dbeafe;
    }

    .catalog-row-extra td {
        background-color: #fee2e2;
        color: #7f1d1d;
    }

    .catalog-row-extra .catalog-rule-id {
        background-color: #fee2e2;
    }

    .catalog-diff-cell {
        background-color: #fef3c7 !important;
        color: #713f12 !important;
    }

    .catalog-version {
        padding: 4px 0;
    }

    .catalog-version-internal {
        border-top: 1px solid #f59e0b;
    }

    .catalog-harmonize {
        margin-top: 4px;
        font-size: 0.75rem;
        font-weight: 700;
        color: #92400e;
    }

    .catalog-row-note {
        margin-top: 6px;
        padding: 5px 7px;
        border: 1px solid currentColor;
        border-radius: 4px;
        font-weight: 700;
    }

    /* Style pour les lignes du tableau */
    .stExpander {
        border: none !important;
        box-shadow: none !important;
    }
    
    .stExpander > div:first-child {
        border-radius: 4px !important;
        background-color: #f5f5f5 !important;
    }
</style>
""", unsafe_allow_html=True)

# Pied de page
st.markdown("---")
st.markdown("Application développée pour l'analyse des formules caryotypiques selon les normes ISCN 2024")
