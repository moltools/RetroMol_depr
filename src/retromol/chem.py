"""
Chemistry module.
"""
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
    def __init__(self, name: str, core: bool, smarts: str) -> None:
        """
        Create a molecular pattern.

        :param str name: Name of molecular pattern.
        :param bool core: Whether or not molecular pattern is a core unit.
        :param str smarts: SMARTS string of molecular pattern.
        :returns: None
        """
        self.name = name 
        self.core = core
        self.smarts = smarts 
        self.compiled = self._compile_pattern(self.smarts)

    def is_core(self) -> bool:
        """
        Check if molecular pattern is a core unit.
        
        :returns: Whether or not molecular pattern is a core unit.
        :rtype: bool
        """
        return self.core

    def _compile_pattern(self, smarts: str) -> Chem.Mol:
        """
        Compile a molecular pattern.
        
        :param str smarts: SMARTS string of molecular pattern.
        :returns: Compiled molecular pattern.
        :rtype: Chem.Mol
        """
        return Chem.MolFromSmarts(smarts)
    
class ReactionRule: 
    def __init__(self, name: str, smirks: str) -> None:
        """
        Create a reaction rule.
        
        :param str name: Name of reaction rule.
        :param str smirks: SMIRKS string of reaction rule.
        :returns: None
        """
        self.name = name 
        self.backward_smirks = smirks 
        self.forward_smirks = self._reverse_smirks(self.backward_smirks)
        self.backward_compiled = self._compile_rule(self.backward_smirks)
        self.forward_compiled = self._compile_rule(self.forward_smirks)

    def _reverse_smirks(self, smirks: str) -> str:
        """
        Reverse a SMIRKS string.
        
        :param str smirks: SMIRKS string.
        :returns: Reversed SMIRKS string.
        :rtype: str
        """
        return ">>".join(reversed(smirks.split(">>")))

    def _compile_rule(self, smirks: str) -> Reaction:
        """
        Compile a reaction rule.
        
        :param str smirks: SMIRKS string of reaction rule.
        :returns: Compiled reaction rule.
        :rtype: Reaction
        """
        return ReactionFromSmarts(smirks)
    
class Molecule:
    def __init__(self, name: str, smiles: str) -> None:
        """
        Create a molecule.
        
        :param str name: Name of molecule.
        :param str smiles: SMILES string of molecule.
        :returns: None
        """
        self.name = name 
        self.smiles = smiles 
        self.compiled = self._compile_molecule(self.smiles)

    def _compile_molecule(self, smiles: str) -> Chem.Mol:
        """
        Compile a molecule.
        
        :param str smiles: SMILES string of molecule.
        :returns: Compiled molecule.
        :rtype: Chem.Mol
        """
        return Chem.MolFromSmiles(smiles)
    
    def apply_rules(
        self,
        reactions: ty.List[ReactionRule]
    ) -> ty.Tuple[Chem.Mol, Tree, ReactionTreeMapping]:
        """
        Apply reaction rules to a molecule.

        :param ty.List[ReactionRule] reactions: List of reaction rules.
        :returns: Reaction products, reaction tree, and mapping of reaction products to molecules.
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

                for results in reaction.backward_compiled.RunReactants([current]):
                    if len(results) > 0:
                        
                        reaction_products = list()

                        for result in results:
                            try:
                                Chem.SanitizeMol(result)
                                # TODO: log sanitization errors.
                            except:
                                # Unable to sanitize molecule. Usually happens because of
                                # explicit valences not being correct, when reaction rule
                                # had unintended side effect on molecule.
                                continue

                            result_encoding = mol_to_encoding(result, N, radius, num_bits)
                            reaction_products.append(result_encoding)

                            if result_encoding not in tree:
                                mols.append(result)

                        tree[current_encoding][reaction.name].add(frozenset(reaction_products))
        
        return self.compiled, tree, mapping

def mol_to_fingerprint(mol: Chem.Mol, radius: int, num_bits: int) -> np.array:
    """
    Turn a molecule into a fingerprint.

    :param Chem.Mol mol: Molecule.
    :param int radius: Fingerprint radius.
    :param int num_bits: Number of bits in fingerprint.
    :returns: Fingerprint.
    :rtype: np.array
    """
    fp_arr = np.zeros((0,), dtype=np.int8)
    fp_vec = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=num_bits)

    DataStructs.ConvertToNumpyArray(fp_vec, fp_arr)

    return fp_arr

def tanimoto_similarity(fp1: np.array, fp2: np.array) -> float:
    """
    Calculate Tanimoto similarity between two fingerprints.

    :param np.array fp1: First fingerprint.
    :param np.array fp2: Second fingerprint.
    :returns: Tanimoto similarity.
    :rtype: float
    """
    return np.logical_and(fp1, fp2).sum() / np.logical_or(fp1, fp2).sum()

def mol_to_encoding(
    mol: Chem.Mol, 
    N: int, 
    radius: int, 
    num_bits: int
) -> np.array:
    """
    Create int hash of reaction product of molecule.

    :param Chem.Mol mol: Molecule.
    :param int N: Number of atoms in molecule.
    :param int radius: Fingerprint radius.
    :param int num_bits: Number of bits in fingerprint.
    :returns: Int hash of reaction product of molecule.
    :rtype: int
    """
    amns = [
        atom.GetIsotope() 
        for atom in mol.GetAtoms() 
        if atom.GetIsotope() > 0
    ]
    amns = np.array([1 if x in amns else 0 for x in np.arange(N)])
    fp = mol_to_fingerprint(mol, radius, num_bits)

    return hash(np.hstack([fp, amns]).data.tobytes())

def identify_mol(
    mol: Chem.Mol, 
    monomers: ty.List[MolecularPattern]
) -> ty.Optional[str]:
    """
    Identify molecule.
    
    :param Chem.Mol mol: Molecule.
    :param ty.List[MolecularPattern] monomers: List of molecular patterns.
    :returns: Name of identified monomer.
    :rtype: ty.Optional[str]
    """
    for monomer in monomers:
        if mol.HasSubstructMatch(monomer.compiled):
            return monomer.name

    return None

def greedy_max_set_cover(
    mol: Chem.Mol, 
    identified: ty.List[ty.Tuple[int, str]],
    mapping: ReactionTreeMapping
) -> ty.List[ty.Tuple[int, str]]:
    """
    Return a greedy maximum set cover of identified monomers.

    :param Chem.Mol mol: Molecule.
    :param ty.List[ty.Tuple[int, str]] identified: List of identified monomers.
    :param ReactionTreeMapping mapping: Mapping of reaction products to molecules.
    :returns: List of selected monomers.
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