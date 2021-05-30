import mols2grid
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt

st.title("Series4 new proposed compounds")

@st.cache(allow_output_mutation=True)
def download_dataset():
    """Loads once then cached for subsequent runs"""
    df = pd.read_csv("data/ReinventResults/processed_results/all_batches.csv")
    df=df.round({"total_score":3})
    return df

# Copy the dataset so any changes are not applied to the original cached version
df = download_dataset().copy()


Score = st.sidebar.slider(
    label="Show compounds with score over",
    min_value=0.,
    max_value=1.,
    value=0.6,
    step=0.01,
)

Activity = st.sidebar.slider(
    label="Show compounds with activity over",
    min_value=0.,
    max_value=1.,
    value=0.6,
    step=0.01,
)


MW = st.sidebar.slider(
    label="Show compounds with molecular weight over",
    min_value=200.,
    max_value=800.,
    value=0.6,
    step=0.01,
)

df_result = df[(df["total_score"] > Score) & (df["MW"] > MW) & (df["ActivityRegression"] > Activity)]
st.write("### EOS-#batch-#compound", df_result)


raw_html = mols2grid.display(df_result.head(100), subset=["eosID","img", "total_score"], tooltip=["eosID", "SMILES","InchiKey", "MW"], selection=False, n_cols=4)._repr_html_()
components.html(raw_html, width=900, height=900, scrolling=True)