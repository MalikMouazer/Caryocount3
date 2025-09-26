import streamlit as st
import pandas as pd
import re
import base64
import io
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
        anom = row['Anomalie']
        type_anom = row['Type']
        score = row['Score ISCN 2024']
        clones = row['Clones']
        explication = row['Explication']
        
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
    """Renvoie un HTML condensé pour la liste des anomalies"""
    html = ""
    for _, row in anomalies_df.iterrows():
        score = row['Score ISCN 2024']
        anomalie = row['Anomalie']
        clones_list = row['Clones'].split(', ')
        clones_clean = list(dict.fromkeys(clones_list))
        clones = ', '.join(clones_clean)
        explication = row['Explication']

        if score == 2:
            color = "#FF5733"
            score_text = "2pts"
        elif score == 1:
            color = "#33A1FF"
            score_text = "1pt"
        else:
            color = "#AAAAAA"
            score_text = "0pt"

        html += f"""
        <div style="margin: 2px 0; padding: 4px 8px; border-left: 3px solid {color}; background-color: #f9f9f9; font-size: 14px;">
            <span style="font-weight: bold;">{clones}</span>
            <span style="color: {color}; font-weight: bold;">[{anomalie}]</span>
            <span style="background-color: #555; color: white; border-radius: 8px; padding: 1px 6px; font-size: 12px;">{score_text}</span>
            <span style="color: #666; margin-left: 8px;">{explication}</span>
        </div>
        """
    return html

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

                # Affichage du tableau avec info-bulles
                st.markdown("### Détail des anomalies")

                # Formatage des anomalies pour l'affichage
                anomalies_df = df.iloc[:-1]  # Exclure la ligne TOTAL
                anomalies_html = format_anomalies_compact(anomalies_df)
                st.markdown(anomalies_html, unsafe_allow_html=True)

                # Affichage du total
                st.markdown(f"**Score total ISCN: {totals['iscn']}**")
                st.markdown(f"**Score total Jondreville: {totals['jondroville']}**")
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
                    if norm == 'count_i':
                        count_i_col = col
                    elif norm == 'count_j':
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
                        anomalies_detail = ", ".join([
                            f"{row_detail['Anomalie']} ({row_detail['Type']}): {row_detail['Score ISCN 2024']} pts"
                            for _, row_detail in anomalies_df.iterrows()
                        ])

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

                display_config = [
                    ("Ligne", 1),
                    ("Formule", 3),
                    ("Comptage ISCN", 1),
                ]
                if has_count_i:
                    display_config.append(("Ref ISCN", 1))
                display_config.append(("Comptage Jon", 1))
                if has_count_j:
                    display_config.append(("Ref Jon", 1))
                display_config.append(("Anomalies détectées", 4))

                # Affichage des résultats
                st.markdown("### Résultats de l'analyse")

                header_cols = st.columns([cfg[1] for cfg in display_config])
                for col_obj, (label, _) in zip(header_cols, display_config):
                    with col_obj:
                        st.markdown(f"**{label}**")

                # Afficher chaque ligne
                for i, (_, row_data) in enumerate(results_df.iterrows()):
                    anomalies = all_anomalies_details[i]
                    matches = match_details[i]

                    line_cols = st.columns([cfg[1] for cfg in display_config])
                    for col_obj, (label, _) in zip(line_cols, display_config):
                        with col_obj:
                            if label == "Comptage ISCN":
                                text = format_display(row_data[label])
                                if matches["iscn"] is not None:
                                    text = f"{text} {'✅' if matches['iscn'] else '❌'}"
                                st.markdown(text)
                            elif label == "Comptage Jon":
                                text = format_display(row_data[label])
                                if matches["jon"] is not None:
                                    text = f"{text} {'✅' if matches['jon'] else '❌'}"
                                st.markdown(text)
                            elif label == "Anomalies détectées":
                                if anomalies["error"]:
                                    st.error(anomalies["message"])
                                else:
                                    html = format_anomalies_compact(anomalies["df"])
                                    st.markdown(html, unsafe_allow_html=True)
                            else:
                                st.markdown(format_display(row_data.get(label)))

                    st.markdown("---")

                # Statistiques de correspondance
                if has_count_i or has_count_j:
                    if has_count_i:
                        total_i = sum(1 for m in match_details if m["iscn"] is not None)
                        match_i = sum(1 for m in match_details if m["iscn"])
                        if total_i:
                            st.success(
                                f"Correspondance ISCN: {match_i}/{total_i} ({int(match_i/total_i*100)}%)"
                            )
                    if has_count_j:
                        total_j = sum(1 for m in match_details if m["jon"] is not None)
                        match_j = sum(1 for m in match_details if m["jon"])
                        if total_j:
                            st.success(
                                f"Correspondance Jondreville: {match_j}/{total_j} ({int(match_j/total_j*100)}%)"
                            )

                # Option d'export Excel
                st.subheader("Exporter les résultats")
                st.markdown(get_excel_download_link(results_df), unsafe_allow_html=True)
                
        except Exception as e:
            st.error(f"Erreur lors de l'analyse du fichier: {str(e)}")

# CSS pour améliorer l'apparence
st.markdown("""
<style>
    .download-button {
        display: inline-block;
        padding: 10px 20px;
        background-color: #4CAF50;
        color: white;
        text-decoration: none;
        border-radius: 4px;
        margin-top: 10px;
        font-weight: bold;
        text-align: center;
    }
    
    .download-button:hover {
        background-color: #45a049;
    }
    
    h3 {
        margin-top: 30px;
        margin-bottom: 20px;
        color: #1E3A8A;
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