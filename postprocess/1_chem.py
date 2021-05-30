# conda activate autocheml
from rdkit import Chem
from rdkit.Chem.Descriptors import qed
from rdkit.Chem.Descriptors import MolLogP
from rdkit.Chem.rdMolDescriptors import CalcNumRings
import pandas as pd
from tqdm import tqdm
from __init__ import FILE
df = pd.read_csv(FILE)
smiles = list(df["Smiles"])

pattern = Chem.MolFromSmarts("*-c1nnc2cncc(-*)n12")
def is_triazo(mol):
    if mol.HasSubstructMatch(pattern):
        return 1
    else:
        return 0

triazo = [is_triazo(Chem.MolFromSmiles(smi)) for smi in tqdm(smiles)]
qeds = [qed(Chem.MolFromSmiles(smi)) for smi in tqdm(smiles)]
slogps = [MolLogP(Chem.MolFromSmiles(smi)) for smi in tqdm(smiles)]
numrings = [CalcNumRings(Chem.MolFromSmiles(smi)) for smi in tqdm(smiles)]

df_ = pd.DataFrame({
    "IsTriazolo": triazo,
    "QED": qeds,
    "SLogP": slogps,
    "NumRings": numrings
})

df_.to_csv("data_1.csv", index=False)
