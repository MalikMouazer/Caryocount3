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
            df, total, error = analyser_formule(formule)
            if error:
                st.error(error)
            else:
                st.success(f"Nombre total d'anomalies détectées: {total}")
                
                # Affichage du tableau avec info-bulles
                st.markdown("### Détail des anomalies")
                
                # Formatage des anomalies pour l'affichage
                anomalies_df = df.iloc[:-1]  # Exclure la ligne TOTAL
                anomalies_html = format_anomalies_compact(anomalies_df)
                st.markdown(anomalies_html, unsafe_allow_html=True)
                
                # Affichage du total
                st.markdown(f"**Score total: {total}**")
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
                # Vérifier si la colonne Count existe
                has_count = 'Count' in df_input.columns
                
                # Création du DataFrame de résultats
                results = []
                all_anomalies_details = []

                for idx, row in df_input.iterrows():
                    formule_fichier = row['Formule']
                    count_manuel = row['Count'] if has_count else None

                    df_analyse, count_auto, error = analyser_formule(formule_fichier)
                    
                    if error:
                        anomalies_detail = error
                        match = "❌" if has_count else "N/A"
                        all_anomalies_details.append({"error": True, "message": error})
                    else:
                        # Extraction des détails des anomalies
                        anomalies_df = df_analyse.iloc[:-1]  # Exclure la ligne TOTAL
                        
                        # Stocker les détails pour l'affichage
                        all_anomalies_details.append({"error": False, "df": anomalies_df})
                        
                        # Texte simple pour l'export
                        anomalies_detail = ", ".join([
                            f"{row['Anomalie']} ({row['Type']}): {row['Score ISCN 2024']} pts"
                            for _, row in anomalies_df.iterrows()
                        ])
                        
                        # Vérification de la correspondance si Count est disponible
                        match = "✅" if has_count and count_auto == count_manuel else "❌" if has_count else "N/A"
                    
                    result_row = {
                        "Ligne": idx + 1,  # +1 pour ne pas inclure l'en-tête du fichier
                        "Formule": formule_fichier,
                        "Comptage automatique": count_auto if not error else "Erreur",
                        "Anomalies détectées": anomalies_detail  # Version texte pour l'export
                    }
                    
                    # Ajouter le comptage manuel si disponible
                    if has_count:
                        result_row["Comptage manuel"] = count_manuel
                        result_row["Correspondance"] = match
                    
                    results.append(result_row)
                
                # Création du DataFrame de résultats
                results_df = pd.DataFrame(results)
                
                # Affichage des résultats
                st.markdown("### Résultats de l'analyse")
                
                # Créer un en-tête de tableau personnalisé
                cols = st.columns([1, 3, 1, 1, 1, 4] if has_count else [1, 3, 1, 6])
                with cols[0]:
                    st.markdown("**Ligne**")
                with cols[1]:
                    st.markdown("**Formule**")
                with cols[2]:
                    st.markdown("**Comptage auto**")
                if has_count:
                    with cols[3]:
                        st.markdown("**Comptage manuel**")
                    with cols[4]:
                        st.markdown("**Correspondance**")
                with cols[-1]:
                    st.markdown("**Anomalies détectées**")
                
                # Afficher chaque ligne avec un expander pour les anomalies
                for i, (_, row_data) in enumerate(results_df.iterrows()):
                    anomalies = all_anomalies_details[i]

                    # Créer une ligne de tableau
                    cols = st.columns([1, 3, 1, 1, 1, 4] if has_count else [1, 3, 1, 6])

                    # Numéro de ligne
                    with cols[0]:
                        st.markdown(f"{row_data['Ligne']}")

                    # Formule
                    with cols[1]:
                        st.markdown(f"{row_data['Formule']}")

                    # Comptage automatique
                    with cols[2]:
                        st.markdown(f"{row_data['Comptage automatique']}")

                    # Comptage manuel et correspondance (si disponible)
                    if has_count:
                        with cols[3]:
                            st.markdown(f"{row_data['Comptage manuel']}")
                        with cols[4]:
                            st.markdown(f"{row_data['Correspondance']}")

                    # Anomalies détectées avec expander
                    with cols[-1]:
                        if anomalies["error"]:
                            st.error(anomalies["message"])
                        else:
                            html = format_anomalies_compact(anomalies["df"])
                            st.markdown(html, unsafe_allow_html=True)
                                        
                    # Ligne de séparation
                    st.markdown("---")
                
                # Statistiques si Count est disponible
                if has_count:
                    nb_total = len(results_df)
                    nb_match = results_df['Correspondance'].value_counts().get("✅", 0)
                    st.success(f"Correspondance: {nb_match}/{nb_total} ({int(nb_match/nb_total*100 if nb_total else 0)}%)")
                
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