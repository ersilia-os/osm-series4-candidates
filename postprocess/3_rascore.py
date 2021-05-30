# conda activate rascore

import pandas as pd
from tqdm import tqdm
from __init__ import FILE

df = pd.read_csv(FILE)
smiles = list(df["Smiles"])

from RAscore import RAscore_NN

scorer = RAscore_NN.RAScorerNN()

scores = [scorer.predict(smi) for smi in tqdm(smiles)]

df_ = pd.DataFrame({
    "RAScore": scores
})

df_.to_csv("data_3.csv", index=False)
