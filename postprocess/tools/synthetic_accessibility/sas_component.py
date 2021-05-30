# conda activate autocheml

import joblib
from typing import List, Tuple
import os

import numpy as np
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem, Descriptors

from .sascorer import calculateScore


PATH = os.path.dirname(os.path.abspath(__file__))

class SASPredict(object):

    def __init__(self):
        self.activity_model = joblib.load(os.path.join(PATH, "SA_score_prediction.pkl"))

    def predict_from_molecules(self, molecules: List) -> np.array:
        if len(molecules) == 0:
            return np.array([])

        descriptors = self._calculate_descriptors(molecules)
        sas_predictions = self.activity_model.predict_proba(descriptors)

        return sas_predictions[:, 1]

    def _calculate_descriptors(self, molecules: List) -> List:
        fingerprints = self._mols_to_fingerprint(molecules)
        descriptors = []

        for idx, mol in enumerate(molecules):
            others = np.array([calculateScore(mol), Descriptors.ExactMolWt(mol)])
            prop_array = np.concatenate([others, fingerprints[idx]]).reshape((1, -1))[0]
            descriptors.append(prop_array)
        return descriptors

    def _mols_to_fingerprint(self, mols) -> List:
        fingerprints = [AllChem.GetHashedMorganFingerprint(mol, 3, nBits=4096) for mol in mols]
        fp_array = []

        for fp in fingerprints:
            numpy_fingreprint = np.zeros((1,))
            DataStructs.ConvertToNumpyArray(fp, numpy_fingreprint)
            fp_array.append(numpy_fingreprint)

        return fp_array

    def _get_props(self, mol):
        molwt = Descriptors.ExactMolWt(mol)

        return molwt

    def _predict_sas(self, smiles: List[str], parameters: dict) -> Tuple[np.array, List]:
        fps, valid_idx = self._smiles_to_fingerprints(smiles, parameters)

        if len(valid_idx) == 0:
            return np.array([]), valid_idx
        activity = self.activity_model.predict_proba(fps, parameters)
        return activity, valid_idx
