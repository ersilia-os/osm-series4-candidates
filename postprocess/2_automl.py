# conda activate autocheml
from tools.autocheml.predict.predict import Predictor
import os

import pandas as pd
from tqdm import tqdm
from __init__ import FILE

df = pd.read_csv(FILE)
smiles = list(df["Smiles"])

def chunker(l, n=10000):
    for i in range(0, len(l), n):
        yield l[i:i + n]

PATH = os.path.abspath("models/autocheml/output_clf")
clf_score = []
for chunk in tqdm(chunker(smiles)):
    clf_score += list(Predictor(chunk, PATH).predict())

PATH = os.path.abspath("models/autocheml/output_reg")
reg_score = []
for chunk in tqdm(chunker(smiles)):
    reg_score += list(Predictor(chunk, PATH).predict())

df_ = pd.DataFrame({
    "ActivityClfAutoML": clf_score,
    "ActivityRegAutoML": reg_score
})

df_.to_csv("data_2.csv", index=False)
