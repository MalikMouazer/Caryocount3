import streamlit as st
import pandas as pd
import re
import base64
import io
import html
import openpyxl
from My_expert_karyo_functions import analyser_formule
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
        f"https://docs.google.com/spreadsheets/d/{sheet_id}/export?format=csv&gid={gid}"
    )

    try:
        df = pd.read_csv(csv_url)
        return df, None
    except Exception as e:
        return None, str(e)


# Configuration de la page
st.set_page_config(
    page_title="Analyseur de Caryotypes",
    page_icon="🧬",
    layout="wide"
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

    def build_pill_text(label, score, explanation, detail_text):
        detail = detail_text
        if explanation != '—':
            detail = f"{detail_text} ({explanation})"
        score_value = render_score_value(score)
        return f"{label} {score_value} : {detail}"

    clone_blocks = {}
    clone_order = []
    has_clone_labels = False

    def add_line(clone_label, line_html):
        key = clone_label or ""
        if key not in clone_blocks:
            clone_blocks[key] = []
            clone_order.append(key)
        clone_blocks[key].append(line_html)

    for _, row in anomalies_df.iterrows():
        score_iscn = row['Score ISCN 2024']
        score_jon = row['Score Jondreville 2020']
        anomalie = html.escape(str(row['Anomalie']))
        type_text = clean_text(row.get('Type'))
        if type_text == '—':
            type_text = anomalie

        explication_iscn = clean_text(row.get('Explication'))
        explication_jon = clean_text(row.get('Explication Jondreville 2020'))
        clone_details = row.get('CloneDetails')
        if not isinstance(clone_details, list) or not clone_details:
            clone_details = [{"label": None, "is_reference": True, "reason": ""}]

        for detail in clone_details:
            label = detail.get('label')
            is_reference = detail.get('is_reference', True)
            reason_raw = detail.get('reason') or ''
            reason_html = clean_text(reason_raw) if reason_raw else ''

            if label:
                has_clone_labels = True

            if is_reference:
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
            line_html = (
                f'<div class="anomaly-line" style="border-left-color: {color};">'
                f'<span class="anomaly-label" style="color: {color};">[{anomalie}]</span>'
                f'<span class="score-pill score-pill-iscn">{build_pill_text("ISCN", line_score_iscn, line_exp_iscn, type_text)}</span>'
                f'<span class="score-pill score-pill-jon">{build_pill_text("Jon", line_score_jon, line_exp_jon, type_text)}</span>'
                '</div>'
            )

            add_line(label, line_html)

    if not clone_blocks:
        return ""

    blocks = []
    show_labels = has_clone_labels and len([lbl for lbl in clone_order if lbl]) >= 2

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


def render_score_totals(score_iscn, score_jon):
    """Crée le bloc HTML présentant les scores globaux."""
    return f"""
    <div class="score-summary-group">
        <div class="score-summary">
            <span class="score-label">ISCN</span>
            <span class="score-pill score-pill-iscn">{render_score_value(score_iscn)}</span>
        </div>
        <div class="score-summary">
            <span class="score-label">Jondreville</span>
            <span class="score-pill score-pill-jon">{render_score_value(score_jon)}</span>
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

# Création des onglets
tab1, tab2 = st.tabs(["Analyse d'une formule", "Analyse d'un fichier"])

# Onglet 1: Analyse d'une formule
with tab1:
    st.subheader("Entrez une formule caryotypique")
    formule = st.text_input("Formule ISCN", placeholder="Ex: 47,XX,+8[20]")
    
    if st.button("Analyser la formule", key="analyser_formule"):
        if formule:
            df, totals, error = analyser_formule(formule)
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
    test_button = st.button("Analyser le fichier de tests")

    uploaded_file = st.file_uploader(
            "Choisir un fichier CSV ou Excel", type=["csv", "xlsx", "xls"], key="file"
        )
 

    df_input = None
    if test_button:
        TEST_SHEET_URL = (
            "https://docs.google.com/spreadsheets/d/"
            "1_QPAzEr3PNHaNVu8Qxuv_hWCUoBQqqRcNSnX-Q0Egm8/edit?gid=98233212#gid=98233212"
        )
        df_input, err = load_google_sheet(TEST_SHEET_URL)
        if err:
            st.error(f"Erreur lors du chargement du fichier de tests : {err}")
            st.stop()
    elif uploaded_file is not None:
        # Déterminer le type de fichier
        if uploaded_file.name.endswith('.csv'):
            df_input = pd.read_csv(uploaded_file)
        else:  # Excel
            df_input = pd.read_excel(uploaded_file)



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

                    df_analyse, totals, error = analyser_formule(formule_fichier)

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
                                clone_details = [{"label": None, "is_reference": True, "reason": ""}]

                            multiple_clones = len([cd for cd in clone_details if cd.get('label')]) >= 2

                            for detail in clone_details:
                                is_reference = detail.get('is_reference', True)
                                reason = clean_detail(detail.get('reason')) if detail.get('reason') else ''
                                clone_label = detail.get('label') if multiple_clones else ''

                                line_score_iscn = row_detail['Score ISCN 2024'] if is_reference else 0
                                line_exp_iscn = exp_iscn if is_reference else (reason or exp_iscn)
                                line_score_jon = row_detail['Score Jondreville 2020'] if is_reference else 0
                                line_exp_jon = exp_jon if is_reference else (reason or exp_jon)

                                prefix = f"{clone_label}: " if clone_label else ""
                                chunk = (
                                    f"{prefix}{type_label}: ISCN {line_score_iscn}"
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
                st.markdown("### Résultats de l'analyse")

                header_html = "".join(
                    f"<th>{html.escape(label)}</th>" for label in display_labels
                )

                body_rows = []
                for i, (_, row_data) in enumerate(results_df.iterrows()):
                    anomalies = all_anomalies_details[i]
                    matches = match_details[i]
                    cells = []

                    for label in display_labels:
                        if label == "Comptage ISCN":
                            raw_value = row_data.get(label)
                            icon = ""
                            if matches["iscn"] is not None:
                                icon = f"<span class='pill-icon'>{'✅' if matches['iscn'] else '❌'}</span>"
                            badge = f"<span class='score-pill score-pill-iscn'>{render_score_value(raw_value)}{icon}</span>"
                            cells.append(f"<td class='score-cell'>{badge}</td>")
                        elif label == "Comptage Jon":
                            raw_value = row_data.get(label)
                            icon = ""
                            if matches["jon"] is not None:
                                icon = f"<span class='pill-icon'>{'✅' if matches['jon'] else '❌'}</span>"
                            badge = f"<span class='score-pill score-pill-jon'>{render_score_value(raw_value)}{icon}</span>"
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
                        else:
                            value = format_display(row_data.get(label))
                            cells.append(f"<td>{html.escape(value)}</td>")

                    body_rows.append(f"<tr>{''.join(cells)}</tr>")

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

                # Statistiques de correspondance + export sur une seule ligne
                columns = st.columns(3)
                renderers = []

                if has_count_i:
                    total_i = sum(1 for m in match_details if m["iscn"] is not None)
                    match_i = sum(1 for m in match_details if m["iscn"])
                    if total_i:
                        msg_i = f"Correspondance ISCN: {match_i}/{total_i} ({int(match_i/total_i*100)}%)"

                        def render_iscn(col, message=msg_i):
                            col.success(message)

                        renderers.append(render_iscn)

                if has_count_j:
                    total_j = sum(1 for m in match_details if m["jon"] is not None)
                    match_j = sum(1 for m in match_details if m["jon"])
                    if total_j:
                        msg_j = f"Correspondance Jondreville: {match_j}/{total_j} ({int(match_j/total_j*100)}%)"

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
        max-height: 500px;
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
