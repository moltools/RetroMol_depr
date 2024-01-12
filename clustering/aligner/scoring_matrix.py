from __future__ import annotations
from enum import Enum, auto
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np

from rdkit import Chem
from rdkit.Chem import rdmolops, AllChem, DataStructs

SCORES = """
A,B,C,D,AA_1,AA_2,AA_3,AA_4,AA_5,AA_6,AA_7,GAP
A,2,-1,-1,-1,-2,-2,-2,-2,-2,-2,-2,-4
B,-1,2,-1,-1,-2,-2,-2,-2,-2,-2,-2,-4
C,-1,-1,2,-1,-2,-2,-2,-2,-2,-2,-2,-4
D,-1,-1,-1,2,-2,-2,-2,-2,-2,-2,-2,-4
AA_1,-2,-2,-2,-2,2,-1,-1,-1,-1,-1,-1,-4
AA_2,-2,-2,-2,-2,-1,2,-1,-1,-1,-1,-1,-4
AA_3,-2,-2,-2,-2,-1,-1,2,-1,-1,-1,-1,-4
AA_4,-2,-2,-2,-2,-1,-1,-1,2,-1,-1,-1,-4
AA_5,-2,-2,-2,-2,-1,-1,-1,-1,2,-1,-1,-4
AA_6,-2,-2,-2,-2,-1,-1,-1,-1,-1,2,-1,-4
AA_7,-2,-2,-2,-2,-1,-1,-1,-1,-1,-1,2,-4
GAP,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,2
"""

# def mol_to_fingerprint(mol: Chem.Mol, n: int = 2048) -> np.array:
#     fingerprint = np.zeros((0,), dtype=int)
#     DataStructs.ConvertToNumpyArray(
#         AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=n), 
#         fingerprint
#     )
#     return fingerprint

# def tanimoto_similarity(fp1: np.array, fp2: np.array) -> float:
#     return (np.logical_and(fp1, fp2).sum() / (float(np.logical_or(fp1, fp2).sum())))

def parse_scoring_matrix() -> Dict[PKModule, Dict[PKModule, int]]:
    pairwise_scores = {}
    scores_str = SCORES
    lines = scores_str.strip().split("\n")
    header = lines[0].split(",")
    for idx, line in enumerate(lines[1:]):
        k_a = PKModule[header[idx]]
        pairwise_scores[k_a] = {}
        k_b = [PKModule[pk] for pk in header]
        # Parse scores line -- first item is PK module name:
        scores = map(int, line.strip().split(",")[1:])
        for k, score in zip(k_b, scores):
            pairwise_scores[k_a][k] = score
    return pairwise_scores


class PKModule(Enum):
    A = auto()
    B = auto()
    C = auto()
    D = auto()
    AA_1 = auto()
    AA_2 = auto()
    AA_3 = auto()
    AA_4 = auto()
    AA_5 = auto()
    AA_6 = auto()
    AA_7 = auto()
    GAP = auto()  # Denotes gap in PK module sequence

    def display_name(self):
        if self == PKModule.GAP:
            return '-'
        else:
            return self.name

    def display_in_alignment(self):
        max_length = max([len(m.display_name()) for m in PKModule])
        padding = max_length - len(self.display_name())
        return (' ' * padding) + self.display_name()

    def logo_color(self):
        cm = plt.get_cmap('gist_rainbow')
        colors = {
            module.display_name(): cm(1. * module_idx/len(PKModule))
            for module_idx, module in enumerate(PKModule)
        }
        return colors[self.display_name()]


class ScoringMatrix:
    def __init__(self) -> None:
        self.pairwise_scores = parse_scoring_matrix()

    # def __call__(self, module1: PKModule, module2: PKModule) -> float:
    #     try:  # Validate scoring matrix on first call
    #         if self.__class__.__call__._called:
    #             return self._score(module1, module2)  # Score modules
    #     except AttributeError:
    #         if not self._validate():  # Validate scoring matrix
    #             raise ValueError('ScoringMatrix is invalid')
    #         self.__class__.__call__._called = True
    #         # Call scoring matrix again after validation
    #         return self(module1, module2)

    # def _validate(self) -> bool:
    #     # Validate if all modules are pairwise represented in scoring matrix
    #      return all([
    #         all([
    #             m2 in self.pairwise_scores[m1]
    #             for m2 in PKModule
    #         ])
    #         for m1 in PKModule
    #     ])

    def score(self, module1: PKModule, module2: PKModule) -> float:
        # try: 
        score = self.pairwise_scores[module1][module2]
        # except KeyError: 
        #     if not isinstance(module1, PKModule) and not isinstance(module2, PKModule):
        #         fp1 = mol_to_fingerprint(Chem.MolFromSmiles(module1))
        #         fp2 = mol_to_fingerprint(Chem.MolFromSmiles(module2))
        #         return (tanimoto_similarity(fp1, fp2) - 0.5) * 20
        #     elif isinstance(module1, PKModule) or isinstance(module2, PKModule):
        #         score = -5.0
        return score