# conda activate autocheml

from tqdm import tqdm
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.cluster import MiniBatchKMeans
from sklearn.neighbors import NearestNeighbors
import pickle
from __init__ import FILE

size = 2048
radius = 3

MAX_N = 200000

df = pd.read_csv(FILE)
smiles = list(df["Smiles"])

useCounts = True
useFeatures = True

def get_fingerprints(smiles):
    fps = [AllChem.GetMorganFingerprint(Chem.MolFromSmiles(mol), radius, useCounts=useCounts, useFeatures=useFeatures) for mol in tqdm(smiles)]
    V = np.zeros((len(fps), size), dtype=np.int32)
    for i, fp in tqdm(enumerate(fps)):
        for idx, v in fp.GetNonzeroElements().items():
            nidx = idx % size
            V[i, nidx] += int(v)
    return V

V = get_fingerprints(smiles)

idxs = np.random.choice([i for i in range(V.shape[0])], min(MAX_N, V.shape[0]), replace=False)

cluster = MiniBatchKMeans(n_clusters=100, verbose=2)
cluster.fit(V[idxs])
labels_100 = cluster.predict(V)

cluster = MiniBatchKMeans(n_clusters=1000, verbose=2)
cluster.fit(V[idxs])
labels_1000 = cluster.predict(V)

smiles_ex = list(pd.read_csv("tools/series4_allsmiles.csv")["canonical"])
V_ex = get_fingerprints(smiles_ex)
neigh = NearestNeighbors(n_neighbors=1, metric="jaccard")
neigh.fit(V_ex)
D, I = neigh.kneighbors(V)
tanimotos = 1 - D[:,0]

df_ = pd.DataFrame(
{
    "Cluster100": labels_100,
    "Cluster1000": labels_1000,
    "TanimotoExisting": tanimotos
}
)

df_.to_csv("data_6.csv", index=False)
