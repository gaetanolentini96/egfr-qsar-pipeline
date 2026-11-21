import streamlit as st
import joblib
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

# Configurazione della pagina Streamlit
st.set_page_config(
    page_title="EGFR QSAR Predictor",
    page_icon="🧬",
    layout="centered"
)

# Caricamento del modello (usiamo @st.cache_resource per caricarlo una sola volta in memoria)
@st.cache_resource
def load_model():
    return joblib.load("egfr_xgboost_model.pkl")

try:
    model = load_model()
except Exception as e:
    st.error(f"Errore nel caricamento del modello: {e}")

# Interfaccia Grafica
st.title("🧬 EGFR (CHEMBL203) Bioactivity Predictor")
st.markdown("""
A production-grade Cheminformatics web app powered by **RDKit** and **XGBoost** 
to predict small-molecule bioactivity against the Epidermal Growth Factor Receptor (EGFR).
""")

st.divider()

# Input per la stringa SMILES (con un esempio di default: Gefitinib)
default_smiles = "COc1cc2c(ncnc2cc1OCCCN3CCOCC3)Nc4ccc(cc4Cl)Cl"
smiles_input = st.text_input(
    "Inserisci la struttura chimica in formato SMILES:", 
    value=default_smiles,
    help="Esempio: Gefitinib (Inibitore EGFR approvato)"
)

if st.button("Esegui Predizione 🚀", type="primary"):
    if not smiles_input:
        st.warning("Per favore inserisci una stringa SMILES valida.")
    else:
        with st.spinner("Analisi molecolare e featurizzazione in corso..."):
            try:
                # 1. Validazione e conversione con RDKit
                mol = Chem.MolFromSmiles(smiles_input)
                if mol is None:
                    st.error("❌ Errore: La stringa SMILES inserita non è valida o non è stata riconosciuta da RDKit.")
                else:
                    # 2. Generazione dei Morgan Fingerprints (Radius=2, 2048 bits)
                    fp = np.array(AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048))
                    
                    # 3. Predizione con il modello XGBoost
                    prediction = model.predict([fp])[0]
                    probability = model.predict_proba([fp])[0][1]
                    
                    # 4. Visualizzazione della molecola
                    img = Draw.MolToImage(mol, size=(500, 350))
                    
                    col1, col2 = st.columns([1, 1])
                    
                    with col1:
                        st.image(img, caption="Struttura Molecolare (RDKit)", use_container_width=True)
                        
                    with col2:
                        st.subheader("Risultato dell'Analisi")
                        if prediction == 1:
                            st.success(f"🟢 **PREDIZIONE: ATTIVO**\n\nIl modello indica un'alta probabilità di inattività/attività biologica mirata ($\le 1000 \\text{{ nM}}$).")
                        else:
                            st.warning(f"🔴 **PREDIZIONE: INATTIVO**\n\nIl composto supera la soglia di attività stabilita ($> 1000 \\text{{ nM}}$).")
                        
                        st.metric(
                            label="Probabilità di Attività (Active Class)", 
                            value=f"{probability:.2%}"
                        )
                        
            except Exception as e:
                st.error(f"Si è verificato un errore durante l'elaborazione: {e}")

st.divider()
st.markdown("*Sviluppato da Gaetano Lentini | Progetto `egfr-qsar-pipeline`*")