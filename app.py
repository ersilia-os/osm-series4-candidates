import mols2grid
import pandas as pd
import streamlit as st
import numpy as np
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt
import matplotlib.pyplot as plt


st.set_page_config(layout="wide")# force wide display
st.title("New candidates for OSM Series 4 by Ersilia Open Source Initiative")

@st.cache
def load_data():
    """Loads once then cached for subsequent runs"""
    df = pd.read_csv("postprocess/210530_EOSI_OSM_Series4_1000.csv")
    return df

df=load_data()

#reorder columns to show activity values first
df=df[["EosId", "ActivityAvg", "MolWeight", "TanimotoExisting", "QED", "SLogP", "SAScore", "RAScore", "IsTriazolo", "NumRings", "ActivityClfAutoML", "ActivityRegAutoML", "ActivityClfGraph", "ActivityRegGraph", "Batch", "Smiles", "InChIKey", "Cluster1000", "Cluster100", ]]

df=df.rename(columns={"Smiles":"SMILES"}) #rename SMILES column so it is recognized by Mols2Vec
df=df.rename(columns={"TanimotoExisting":"Tanimoto"}) #rename SMILES column so it is recognized by Mols2Vec
df=df.round(4) #round activity to 4 decimals


# Create Columns with selection options on the left and the dataframe displayed on the right
col1, col2=st.beta_columns((1,3)) #create two columns, for slider optins (col1) and table (col2)

col1.header("Search by:")
EosId=df["EosId"].tolist()
molecule = col1.selectbox("Molecule",([None] + EosId))

Activity = col1.slider(
    label="Predicted activity against P. falciparum",
    min_value=0.,
    max_value=1.,
    value=0.6,
    step=0.01,
)


MW = col1.slider(
    label="Molecular weight",
    min_value=200.,
    max_value=600.,
    value=300.,
    step=0.01,
)

SA = col1.slider(
    label="Synthetic accessibility",
    min_value=0.,
    max_value=1.,
    value=0.6,
    step=0.01,
)


col2.header("EOS - #batch - #compound") 
if molecule == None:
    df_result = df[(df["MolWeight"] > MW) & (df["ActivityAvg"] > Activity) & (df["SAScore"] > SA)]
    col2.dataframe(df_result,width=None, height=400)
else:
    col2.write(df[df["EosId"] == molecule])

#Display molecular structure of compounds selected above
st.subheader("Molecular Structures")
if molecule == None:
    raw_html = mols2grid.display(df_result.head(100), subset=["EosId","img", "ActivityAvg"], tooltip=["EosId", "SMILES","InChIKey", "Cluster1000"], selection=False, n_cols=6)._repr_html_()
    components.html(raw_html,height=600, scrolling=True)
else:
    raw_html = mols2grid.display(df[df["EosId"] == molecule], subset=["EosId","img", "ActivityAvg"], tooltip=["EosId", "SMILES","InChIKey", "Cluster1000"], selection=False, n_cols=6)._repr_html_()
    components.html(raw_html,height=600, scrolling=True)