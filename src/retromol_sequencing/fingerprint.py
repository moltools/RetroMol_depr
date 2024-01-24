"""
Module for biosynthetic fingerprinting.
"""
import math 
import typing as ty
from enum import Enum

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, GraphDescriptors

class CompoundClassMapping(Enum):
    """
    Enum class that maps motif class to coordinates on unit circle.

    For more information on the compound classes, see:
    https://srobinson.shinyapps.io/AdenylPred/#section-substrate-groups 
    """
    A =                     ( 1.000,  0.000)
    B =                     ( 0.841,  0.541)
    C =                     ( 0.415,  0.910)
    D =                     (-0.142,  0.990)
    PolarAndCharged =       (-0.655,  0.756)
    SmallHydrophobic =      (-0.959,  0.282)
    SmallNonHydrophobic =   (-0.655, -0.756)
    Tiny =                  ( 0.841, -0.541)
    Bulky =                 (-0.959, -0.282)
    CyclicAliphatic =       (-0.142, -0.990)
    SulfurContaining =      ( 0.415, -0.910)
    Undefined =             ( 0.000,  0.000)

    def get_mapping(item: str) -> ty.Tuple[float, float]:
        """
        Get mapping on unit circle for motif class.
        
        :param str item: Item.
        :return: Mapping.
        :rtype: ty.Tuple[float, float]
        :raises ValueError: If item is not a valid motif class.
        """
        return CompoundClassMapping[item].value
    
    @classmethod
    def get_scoring_matrix(cls) -> str:
        """
        Get scoring matrix for motif classes.
        
        :return: Scoring matrix.
        :rtype: str
        """
        return """A,B,C,D,PolarAndCharged,SmallHydrophobic,SmallNonHydrophobic,Tiny,Bulky,CyclicAliphatic,SulfurContaining,Undefined
A,3,2,2,2,0,0,0,0,0,0,0,-2
B,2,3,2,2,0,0,0,0,0,0,0,-2
C,2,2,3,2,0,0,0,0,0,0,0,-2
D,2,2,2,3,0,0,0,0,0,0,0,-2
PolarAndCharged,0,0,0,0,3,2,2,2,2,2,2,-2
SmallHydrophobic,0,0,0,0,2,3,2,2,2,2,2,-2
SmallNonHydrophobic,0,0,0,0,2,2,3,2,2,2,2,-2
Tiny,0,0,0,0,2,2,2,3,2,2,2,-2
Bulky,0,0,0,0,2,2,2,2,3,2,2,-2
CyclicAliphatic,0,0,0,0,2,2,2,2,2,3,2,-2
SulfurContaining,0,0,0,0,2,2,2,2,2,2,3,-2
Undefined,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,3"""
    
def get_biosynthetic_fingerprint(seq: ty.List[str]) -> np.ndarray:
    """
    Create fingerprint from sequence, based on Chaos Game Representation.

    :param ty.List[str] seq: List of motif classes.
    :return: Biosynthetic fingerprint (N=3600).
    :rtype: np.ndarray
    """
    # Make empty 3600 count array.
    fingerprint = np.zeros(3600, dtype=np.int8)

    def find_bin(x: float, y: float) -> int:
        """
        Find bin index for given coordinates.

        :param float x: x-coordinate.
        :param float y: y-coordinate.
        :return: Bin index.
        :rtype: int
        """
        # Get clockwise angle from x-axis.
        angle = math.degrees(math.atan2(y, x))
        if angle < 0: angle += 360
        angle = math.ceil(angle)

        # Get distance to origin and multiply by 10 do get distance to 
        # origin as integer between 0 and 10.
        dist = math.ceil(math.sqrt(x**2 + y**2) * 10)

        # Get bin index.
        bin_index = angle * dist

        return bin_index
    
    prev_loc = (0, 0)
    for motif in seq:
        next_loc = CompoundClassMapping.get_mapping(motif)        
        middle_loc = ((prev_loc[0] + next_loc[0]) / 2, (prev_loc[1] + next_loc[1]) / 2)
        bin_index = find_bin(*middle_loc)
        fingerprint[bin_index] += 1
        prev_loc = middle_loc

    # Normalize fingerprint.
    fingerprint = fingerprint / np.linalg.norm(fingerprint)

    return fingerprint

def get_amino_acid_fingerprint(mol: Chem.Mol) -> np.ndarray:
    """
    Get fingerprint for molecule.
    
    :param Chem.Mol mol: Molecule.
    :returns: Fingerprint.
    :rtype: np.ndarray
    """
    def feat_num_heavy_atoms(mol):
        return float(sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1]))
    
    def feat_num_sulfur_atoms(mol):
        return float(sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [16]))

    def feat_num_3_rings(mol):
        sssr = Chem.GetSSSR(mol)
        return float(sum(1 for ring in sssr if len(ring) == 3))

    def feat_num_5_rings(mol):
        sssr = Chem.GetSSSR(mol)
        return float(sum(1 for ring in sssr if len(ring) == 5))

    def feat_num_6_rings(mol):
        sssr = Chem.GetSSSR(mol)
        return float(sum(1 for ring in sssr if len(ring) == 6))

    # More information: https://www.rdkit.org/docs/source/rdkit.Chem.Lipinski.html
    def feat_num_hydrogen_bond_acceptors(mol):
        return float(Lipinski.NumHAcceptors(mol))

    # More information: https://www.rdkit.org/docs/source/rdkit.Chem.Lipinski.html
    def feat_num_hydrogen_bond_donors(mol):
        return float(Lipinski.NumHDonors(mol))

    # More information: https://www.rdkit.org/docs/source/rdkit.Chem.GraphDescriptors.html
    def feat_bertzct(mol):
        return float(GraphDescriptors.BertzCT(mol))

    # More information: https://www.rdkit.org/docs/source/rdkit.Chem.Lipinski.html
    def feat_num_nhs_ohs(mol):
        return float(Lipinski.NHOHCount(mol))
    
    def feat_mol_weight(mol):
        return float(Descriptors.ExactMolWt(mol))
    
    def feat_logp(mol):
        return float(Crippen.MolLogP(mol))

    return np.array([
        feat_num_heavy_atoms(mol),
        feat_num_sulfur_atoms(mol),
        feat_num_3_rings(mol),
        feat_num_5_rings(mol),
        feat_num_6_rings(mol),
        feat_num_hydrogen_bond_acceptors(mol),
        feat_num_hydrogen_bond_donors(mol),
        feat_bertzct(mol),
        feat_num_nhs_ohs(mol),
        feat_mol_weight(mol),
        feat_logp(mol),
    ])

def amino_acid_class_to_label(classification: int) -> str:
    """
    Convert amino acid class to label.

    :param int classification: Classification.
    :return: Label.
    :rtype: str
    """
    return {
        0: "PolarAndCharged",
        1: "SmallHydrophobic",
        2: "Bulky",
        3: "SmallNonHydrophobic",
        4: "CyclicAliphatic",
        5: "Tiny",
        6: "SulfurContaining"
    }[classification]