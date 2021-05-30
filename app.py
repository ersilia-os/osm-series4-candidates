import mols2grid
import pandas as pd
import streamlit as st
import numpy as np
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt
st.set_page_config(layout="wide")# force wide display
st.title("Series4 new candidates")

@st.cache
def load_data():
    """Loads once then cached for subsequent runs"""
    df = pd.read_csv("data/ReinventResults/ProcessedResults/sample.csv")
    return df

df=load_data()

df=df.round(3) #round numbers to 3 decimals
df=df.rename(columns={"Smiles":"SMILES"}) #rename SMILES column so it is recognized by Mols2Vec


col1, col2=st.beta_columns((1,3)) #create two columns, for slider optins (col1) and table (col2)

col1.header("Search by:")


Activity = col1.slider(
    label="Show compounds with activity over",
    min_value=0.,
    max_value=1.,
    value=0.6,
    step=0.01,
)


MW = col1.slider(
    label="Show compounds with molecular weight over",
    min_value=200.,
    max_value=800.,
    value=300.,
    step=0.01,
)

RA = col1.slider(
    label="Show compounds with RA score over",
    min_value=0.,
    max_value=1.,
    value=0.6,
    step=0.01,
)

col2.header("EOS - #batch - #compound")
            
df_result = df[(df["MolWeight"] > MW) & (df["AvActivity"] > Activity) & (df["RAScore"] > RA)]
col2.dataframe(df_result,width=None, height=600)



st.subheader("Molecular Structures")
raw_html = mols2grid.display(df_result.head(100), subset=["EosId","img", "AvActivity"], tooltip=["EosId", "SMILES","InChIKey", "Cluster1000"], selection=False, n_cols=6)._repr_html_()
components.html(raw_html,height=600, scrolling=True)

st.subheader("Molecular weight distribution")
hist_values=np.histogram(df["MolWeight"])
st.bar_chart(hist_values)