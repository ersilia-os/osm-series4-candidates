# conda activate autocheml

import pandas as pd
import numpy as np
from scipy.stats import rankdata
import collections

df = pd.read_csv("data_7.csv")
df = df.round(3)

columns = [
    "ActivityClfAutoML",
    "ActivityRegAutoML",
    "ActivityClfGraph",
    "ActivityRegGraph"
]

for col in columns:
    values = np.array(df[col])
    ranks = rankdata(values, method="ordinal")/len(values)
    df[col] = ranks

score = []
for r in df[columns].values:
    score += [np.mean(r)]
df["ActivityAvg"] = score

# select 1000 examples
clusters = collections.defaultdict(list)
for i, r in enumerate(df[["Cluster1000", "ActivityAvg"]].values):
    clusters[r[0]] += [(i, r[1])]

idxs = []
for k,v in clusters.items():
    v = sorted(v, key = lambda x: -x[1])
    idxs += [v[0][0]]

df_f = df.iloc[sorted(idxs)]

df.to_csv("210530_EOSI_OSM_Series4_All.csv", index=False)
df_f.to_csv("210530_EOSI_OSM_Series4_1000.csv", index=False)
