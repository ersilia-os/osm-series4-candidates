# conda activate autocheml

import pandas as pd
from __init__ import FILE

data_files = ["data_{0}.csv".format(x) for x in "123456"]

df = pd.read_csv(FILE)

columns = [
    "EosId",
    "Smiles",
    "InChIKey",
    "Batch",
    "IsTriazolo",
    "TanimotoExisting",
    "Cluster100",
    "Cluster1000",
    "MolWeight",
    "SLogP",
    "NumRings",
    "QED",
    "SAScore",
    "RAScore",
    "ActivityClfAutoML",
    "ActivityRegAutoML",
    "ActivityClfGraph",
    "ActivityRegGraph"
]

df = pd.concat([df] + [pd.read_csv(f) for f in data_files], axis=1)

df = df[columns]

df.to_csv("data_7.csv", index=False)
