#!/usr/bin/env python3
"""This module contains the chemistry utilities for the RetroMol package."""

import logging
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


def neutarlize_molecule(mol: Chem.Mol) -> Chem.Mol:
    """Neutralize the molecule.

    :param mol: The molecule.
    :type mol: Chem.Mol
    :return: The neutralized molecule.
    :rtype: Chem.Mol
    """
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol


class MolecularPattern:
    """A class to represent a molecular pattern.

    :param name: The name of the molecular pattern.
    :type name: str
    :param pattern: The SMARTS pattern.
    :type pattern: str
    """

    def __init__(self, name: str, pattern: str) -> None:
        self.name = name
        self.pattern = pattern
        self.compiled = Chem.MolFromSmarts(pattern)


class ReactionRule:
    """A class to represent a reaction rule.

    :param name: The name of the reaction rule.
    :type name: str
    :param pattern: The SMARTS reaction pattern.
    :type pattern: str
    :type properties: ty.Any
    """

    def __init__(self, name: str, pattern: str) -> None:
        self.name = name
        self.patterns = pattern
        self.compiled = ReactionFromSmarts(pattern)


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

    def neutralize(self) -> None:
        """Neutralize the molecule."""
        self.compiled = neutarlize_molecule(self.compiled)

    def apply_rules(
        self, reactions: ty.List[ReactionRule]
    ) -> ty.Tuple[Chem.Mol, Tree, ReactionTreeMapping]:
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

            for reaction in reactions:
                for results in reaction.compiled.RunReactants([current]):
                    if len(results) > 0:
                        reaction_products = list()
                        for result in results:
                            try:
                                Chem.SanitizeMol(result)
                            except Exception:
                                continue

                            result_encoding = mol_to_encoding(
                                result, N, radius, num_bits
                            )
                            reaction_products.append(result_encoding)

                            if result_encoding not in tree:
                                mols.append(result)

                        tree[current_encoding][reaction.name].add(
                            frozenset(reaction_products)
                        )

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


def identify_mol(
    mol: Chem.Mol, monomers: ty.List[MolecularPattern]
) -> ty.Optional[str]:
    """Identify the molecule.

    :param mol: The molecule.
    :type mol: Chem.Mol
    :param monomers: The list of molecular patterns.
    :type monomers: ty.List[MolecularPattern]
    :return: The identifier of the identified molecule.
    :rtype: ty.Optional[str]
    """
    for monomer in monomers:
        if mol.HasSubstructMatch(monomer.compiled):
            return monomer.name

    return None


def greedy_max_set_cover(
    mol: Chem.Mol, identified: ty.List[ty.Tuple[int, str]], mapping: ReactionTreeMapping
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
    logger = logging.getLogger(__name__)
    
    superset = set()
    for atom in mol.GetAtoms():
        if atom.GetIsotope() > 0:
            superset.add(atom.GetIsotope())

    subsets = list()
    for node, node_id in identified:
        submol = mapping[node]
        amns = set(
            [atom.GetIsotope() for atom in submol.GetAtoms() if atom.GetIsotope() > 0]
        )
        subsets.append((amns, (node, node_id)))

    sorted_subsets = sorted(subsets, key=lambda x: len(x[0]), reverse=True)

    logger.debug(f"Founds {len(sorted_subsets)} subsets. Performing greedy maximum set cover ...")

    selected_subsets = []
    covered_elements = set()
    for subset, info in sorted_subsets:
        uncovered_elements = subset - covered_elements

        # Make sure that only a subset is selected if all elements are uncovered.
        if uncovered_elements != subset:
            continue

        if uncovered_elements:
            selected_subsets.append(info)
            covered_elements.update(uncovered_elements)

    logger.debug(f"Performed greedy maximum set cover and selected {len(selected_subsets)} subsets.")

    return selected_subsets
