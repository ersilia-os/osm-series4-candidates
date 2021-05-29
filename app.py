import mols2grid
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt

st.title("Filter new molecules by score")

@st.cache(allow_output_mutation=True)
def download_dataset():
    """Loads once then cached for subsequent runs"""
    df = pd.read_csv("data/processed_results/cleaned1.csv")
    df=df.round({"total_score":3})
    return df

# Copy the dataset so any changes are not applied to the original cached version
df = download_dataset().copy()


score = st.slider(
    label="Show compounds with score over",
    min_value=0.,
    max_value=1.,
    value=0.6,
    step=0.01,
)

df_result = df[df["total_score"] > score]
st.write("### Series 4 derived compounds: eos-#batch-#compound", df_result)


raw_html = mols2grid.display(df_result.head(100), subset=["eosID","img", "total_score"], tooltip=["eosID", "SMILES","InchiKey", "MW"], selection=False, n_cols=7)._repr_html_()
components.html(raw_html, width=1250, height=1250, scrolling=True)