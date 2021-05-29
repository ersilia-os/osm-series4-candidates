import mols2grid
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt

st.title("Filter new molecules by score")

df = pd.read_csv("data/processed_results/cleaned1.csv")
df=df.round({"total_score":3})
score = st.slider(
    label="Show compounds with score over",
    min_value=0.,
    max_value=1.,
    value=0.6,
    step=0.01,
)

df_result = df[df["total_score"] > score]
st.write(df_result)


raw_html = mols2grid.display(df_result, subset=["eosID","img", "total_score"], tooltip=["eosID", "SMILES","InchiKey", "MW"], selection=False)._repr_html_()
components.html(raw_html, width=900, height=900, scrolling=True)