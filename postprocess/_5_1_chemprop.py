from tqdm import tqdm
import pandas as pd
from __init__ import FILE

df = pd.read_csv(FILE)
smiles = list(df["Smiles"])

with open("_chemprop.csv", "w") as f:
    f.write("smiles\n")
    for smi in smiles:
        f.write("{0}\n".format(smi))
