#!/usr/bin/env python3
"""This module contains the chemistry utilities for the RetroMol package."""
import typing as ty
from collections import defaultdict

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.rdChemReactions import ChemicalReaction, ReactionFromSmarts

Tree = ty.Dict[int, ty.Set]
ReactionTreeMapping = ty.Dict[int, Chem.Mol]
MonomerGraphMapping = ty.Dict[int, ty.Tuple[int, str]]
Reaction = ChemicalReaction

class MolecularPattern:
    """A class to represent a molecular pattern.
    
    :param identifier: The identifier of the molecular pattern.
    :type identifier: str
    :param patterns: The list of SMARTS patterns.
    :type patterns: ty.List[str]
    :param properties: The properties of the molecular pattern.
    :type properties: ty.Any
    """
    def __init__(self, identifier: str, patterns: ty.List[str], properties) -> None:
        self.identifier = identifier
        self.patterns = patterns
        self.properties = properties
        self.compiled = [Chem.MolFromSmarts(pattern) for pattern in patterns]

class ReactionRule:
    """A class to represent a reaction rule.
    
    :param identifier: The identifier of the reaction rule.
    :type identifier: str
    :param patterns: The list of SMARTS patterns.
    :type patterns: ty.List[str]
    :param properties: The properties of the reaction rule.
    :type properties: ty.Any
    """
    def __init__(self, identifier: str, patterns: ty.List[str], properties) -> None:
        self.identifier = identifier
        self.patterns = patterns
        self.properties = properties
        self.compiled = [ReactionFromSmarts(pattern) for pattern in patterns]

class Molecule:
    """A class to represent a molecule.
    
    :param name: The name of the molecule.
    :type name: str
    :param smiles: The SMILES representation of the molecule.
    :type smiles: str
    """
    def __init__(self, name: str, smiles: str) -> None:
        self.name = name
        self.smiles = smiles
        self.compiled = Chem.MolFromSmiles(smiles)

    def apply_rules(self, reactions: ty.List[ReactionRule]) -> ty.Tuple[Chem.Mol, Tree, ReactionTreeMapping]:
        """Apply the reaction rules to the molecule.
        
        :param reactions: The list of reaction rules.
        :type reactions: ty.List[ReactionRule]
        :return: The reaction tree.
        :rtype: ty.Tuple[Chem.Mol, Tree, ReactionTreeMapping]
        """
        radius = 2
        num_bits = 2048
        tree = defaultdict(lambda: defaultdict(set))
        mapping = dict()

        N = self.compiled.GetNumAtoms()
        for atom in self.compiled.GetAtoms():
            atom.SetIsotope(atom.GetIdx() + 1)

        mols = [self.compiled]
        while mols:
            current = mols.pop()
            current_encoding = mol_to_encoding(current, N, radius, num_bits)
            tree[current_encoding]
            mapping[current_encoding] = current

            for reaction_rules in reactions:
                for reaction_rule in reaction_rules.compiled:
                    for results in reaction_rule.RunReactants([current]):
                        if len(results) > 0:
                            reaction_products = list()
                            for result in results:
                                try:
                                    Chem.SanitizeMol(result)
                                except Exception:
                                    continue

                                result_encoding = mol_to_encoding(result, N, radius, num_bits)
                                reaction_products.append(result_encoding)

                                if result_encoding not in tree:
                                    mols.append(result)

                            tree[current_encoding][reaction_rules.identifier].add(frozenset(reaction_products))

        return self.compiled, tree, mapping

def mol_to_fingerprint(mol: Chem.Mol, radius: int, num_bits: int) -> np.array:
    """Convert a molecule to a fingerprint.
    
    :param mol: The molecule.
    :type mol: Chem.Mol
    :param radius: The radius of the fingerprint.
    :type radius: int
    :param num_bits: The number of bits.
    :type num_bits: int
    :return: The fingerprint.
    :rtype: np.array
    """
    fp_arr = np.zeros((0,), dtype=np.int8)
    fp_vec = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=num_bits)
    DataStructs.ConvertToNumpyArray(fp_vec, fp_arr)
    return fp_arr

def tanimoto_similarity(fp1: np.array, fp2: np.array) -> float:
    """Calculate the Tanimoto similarity between two fingerprints.
    
    :param fp1: The first fingerprint.
    :type fp1: np.array
    :param fp2: The second fingerprint.
    :type fp2: np.array
    :return: The Tanimoto similarity.
    :rtype: float
    """
    return np.logical_and(fp1, fp2).sum() / np.logical_or(fp1, fp2).sum()

def mol_to_encoding(mol: Chem.Mol, N: int, radius: int, num_bits: int) -> np.array:
    """Convert a molecule to an encoding.

    :param mol: The molecule.
    :type mol: Chem.Mol
    :param N: The number of atoms.
    :type N: int
    :param radius: The radius of the fingerprint.
    :type radius: int
    :param num_bits: The number of bits.
    :type num_bits: int
    :return: The encoding.
    :rtype: np.array
    """
    amns = [atom.GetIsotope() for atom in mol.GetAtoms() if atom.GetIsotope() > 0]
    amns = np.array([1 if x in amns else 0 for x in np.arange(N)])
    fp = mol_to_fingerprint(mol, radius, num_bits)
    return hash(np.hstack([fp, amns]).data.tobytes())

def identify_mol(mol: Chem.Mol, monomers: ty.List[MolecularPattern]) -> ty.Optional[str]:
    """Identify the molecule.
    
    :param mol: The molecule.
    :type mol: Chem.Mol
    :param monomers: The list of molecular patterns.
    :type monomers: ty.List[MolecularPattern]
    :return: The identifier of the identified molecule.
    :rtype: ty.Optional[str]
    """
    for monomer in monomers:
        for compiled_pattern in monomer.compiled:
            if mol.HasSubstructMatch(compiled_pattern):
                return {"identifier": monomer.identifier, "properties": monomer.properties}
    return None

def greedy_max_set_cover(
    mol: Chem.Mol,
    identified: ty.List[ty.Tuple[int, str]],
    mapping: ReactionTreeMapping
) -> ty.List[ty.Tuple[int, str]]:
    """Return a greedy maximum set cover of identified monomers.

    :param mol: The molecule.
    :type mol: Chem.Mol
    :param identified: The list of identified monomers.
    :type identified: ty.List[ty.Tuple[int, str]]
    :param mapping: The mapping of the reaction tree.
    :type mapping: ReactionTreeMapping
    :return: The selected monomers.
    :rtype: ty.List[ty.Tuple[int, str]]
    """
    superset = set()
    for atom in mol.GetAtoms():
        if atom.GetIsotope() > 0:
            superset.add(atom.GetIsotope())

    subsets = list()
    for node, node_id in identified:
        submol = mapping[node]
        amns = set([
            atom.GetIsotope()
            for atom in submol.GetAtoms()
            if atom.GetIsotope() > 0
        ])
        subsets.append((amns, (node, node_id)))

    sorted_subsets = sorted(subsets, key=lambda x: len(x[0]), reverse=True)

    selected_subsets = []
    covered_elements = set()
    for subset, info in sorted_subsets:
        uncovered_elements = subset - covered_elements
        if uncovered_elements:
            selected_subsets.append(info)
            covered_elements.update(uncovered_elements)

    return selected_subsets
