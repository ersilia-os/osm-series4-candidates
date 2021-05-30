# conda activate autocheml

from tools.synthetic_accessibility.sas_component import SASPredict
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from __init__ import FILE

df = pd.read_csv(FILE)
smiles = list(df["Smiles"])

def chunker(l, n=10000):
    for i in range(0, len(l), n):
        yield l[i:i + n]

scorer = SASPredict()
score = []
for chunk in tqdm(chunker(smiles)):
    mols = [Chem.MolFromSmiles(smi) for smi in chunk]
    score += list(scorer.predict_from_molecules(mols))

df_ = pd.DataFrame({
    "SAScore": score,
})

df_.to_csv("data_4.csv", index=False)
