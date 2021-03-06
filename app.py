import mols2grid
import pandas as pd
import streamlit as st
import numpy as np
import streamlit.components.v1 as components


st.set_page_config(layout="wide")# force wide display
st.title("New candidates for OSM Series 4 by Ersilia Open Source Initiative")
st.write("Know about [Ersilia](https://ersilia.io) | Code for [the analysis](https://colab.research.google.com/drive/1MK4UJP6Vw1FjaVaTu9CwpBb2XrjPw9xz?usp=sharing) | Repository of [results](https://github.com/ersilia-os/osm-series4-candidates) | [Tree Map](https://ersilia-os.github.io/osm-series4-candidates/tmap/osm_eosi_s4.html) | Last modified 30th of May 2021")

@st.cache
def load_data():
    """Loads once then cached for subsequent runs"""
    df = pd.read_csv("postprocess/210530_EOSI_OSM_Series4_1000_substr.csv")
    return df

df=load_data()

#reorder columns to show activity values first
df=df[["EosId", "ActivityAvg", "MolWeight", "TanimotoExisting", "QED", "SLogP", "SAScore", "RAScore", "IsTriazolo", "NumRings", "ActivityClfAutoML", "ActivityRegAutoML", "ActivityClfGraph", "ActivityRegGraph", "Batch", "Smiles", "InChIKey", "Cluster1000", "Cluster100","Heteroaryl", "Phenyl", "Para", "Meta", "Orto"]]

df=df.rename(columns={"Smiles":"SMILES"}) #rename SMILES column so it is recognized by Mols2Vec
df=df.rename(columns={"TanimotoExisting":"Tanimoto"}) #rename SMILES column so it is recognized by Mols2Vec
df=df.round(3) #round activity to 4 decimals

#create 5 columns to contain the checkboxes for substructural elements
st.header("Restrict the Triazole substituent to the following substructures")
col1, col2, col3, col4, col5 =st.beta_columns(5) 

check = col1.selectbox('Heteroaryl',["Any", "Yes", "No"], help="Aryl ring with heteroatoms in any of its positions. Does not exclude substituents")
if check == "Any":
    df = df
elif check == "Yes":
    df = df[(df["Heteroaryl"] == 1)]
elif check == "No":
    df = df[(df["Heteroaryl"] == 0)]
    

check2 = col2.selectbox('Phenyl', ["Any", "Yes", "No"], help="Phenyl ring. Does not exclude substituents")
if check2 == "Any":
    df = df
elif check2 == "Yes":
    df = df[(df["Phenyl"] == 1)]
elif check2 == "No":
    df = df[(df["Phenyl"] == 0)] 

check3 = col3.selectbox('Para', ["Any", "Yes", "No"], help="Phenyl ring or heteroaryl with para- substituents. Does not exclude other substituent positions")
if check3 == "Any":
    df = df
elif check3 == "Yes":
    df = df[(df["Para"] == 1)]
elif check3 == "No":
    df = df[(df["Para"] == 0)] 
    
check4 = col4.selectbox('Meta', ["Any", "Yes", "No"], help="Phenyl ring or heteroaryl with meta- substituents. Does not exclude other substituent positions")
if check4 == "Any":
    df = df
elif check4 == "Yes":
    df = df[(df["Meta"] == 1)]
elif check4 == "No":
    df = df[(df["Meta"] == 0)] 
    
check5 = col5.selectbox('Orto', ["Any", "Yes", "No"], help="Phenyl ring or heteroaryl with -orto- substituents. Does not exclude other substituent positions")
if check5 == "Any":
    df = df
elif check5 == "Yes":
    df = df[(df["Orto"] == 1)]
elif check5 == "No":
    df = df[(df["Orto"] == 0)] 

# Create Columns with selection options on the left and the dataframe displayed on the right
col1, col2=st.beta_columns((1,3)) #create two columns, for slider optins (col1) and table (col2)

col1.header("Filter by")
EosId=df["EosId"].tolist()
molecule = col1.selectbox("Molecule",(["All"] + EosId))

Activity = col1.slider(
    label="Predicted activity against P. falciparum",
    min_value=0.,
    max_value=1.,
    value=(0.6, 1.),
    step=0.01,
)

MW = col1.slider(
    label="Molecular weight",
    min_value=200.,
    max_value=600.,
    value=(350.,550.),
    step=0.01,
)

SA = col1.slider(
    label="Synthetic accessibility",
    min_value=0.,
    max_value=1.,
    value=(0.6,1.),
    step=0.01,
)

SL = col1.slider(
    label="Solubility (SLogP)",
    min_value=-1.,
    max_value=10.,
    value=(-1.,4.5),
    step=0.01,
)



col2.header("Selected candidates from pool of 1000 representatives")
if molecule == "All":
    df = df[(df["ActivityAvg"] <= Activity[1]) & (df["ActivityAvg"] >= Activity[0])]
    df = df[(df["MolWeight"] <= MW[1]) & (df["MolWeight"]>= MW[0])]
    df = df[(df["SLogP"] <= SL[1]) & (df["SLogP"]>= SL[0])]
    df_result = df[(df["SAScore"] <= SA[1]) & (df["SAScore"]>= SA[0])]
    col2.dataframe(df_result,width=None, height=500)
else:
    df_result=df[df["EosId"] == molecule]
    col2.write(df_result)

col1.write("{0} molecules selected".format(df_result.shape[0]))

#Display molecular structure of compounds selected above

col1.header("Molecule structures")
n_cols = col2.slider(
    label="Number of molecular structre grid columns to display",
    min_value=4,
    max_value=12,
    value=6,
    step=1)


if molecule == "All":
    raw_html = mols2grid.display(df_result.head(100), subset=["EosId","img", "ActivityAvg"], tooltip=["EosId", "SMILES","InChIKey", "Cluster1000"], selection=False, n_cols=n_cols)._repr_html_()
    components.html(raw_html, width=None, height=600, scrolling=True)
else:
    raw_html = mols2grid.display(df_result, subset=["EosId","img", "ActivityAvg"], tooltip=["EosId", "SMILES","InChIKey", "Cluster1000"], selection=False, n_cols=n_cols)._repr_html_()
    components.html(raw_html,height=600, scrolling=True)

st.write("[Download 1k candidates](https://raw.githubusercontent.com/ersilia-os/osm-series4-candidates/main/postprocess/210530_EOSI_OSM_Series4_1000.csv) | [Download 100k candidates](https://raw.githubusercontent.com/ersilia-os/osm-series4-candidates/main/postprocess/210530_EOSI_OSM_Series4_All.csv) |??Disclaimer: this is a first round of generative models, please give us feedback through OSM GitHub Issue [#34](https://github.com/OpenSourceMalaria/Series4_PredictiveModel/issues/34).")
