import typing as ty 
from collections import defaultdict 
from copy import deepcopy

import numpy as np 
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.rdChemReactions import ChemicalReaction, ReactionFromSmarts 

Tree = ty.Dict[int, ty.Set]
ReactionTreeMapping = ty.Dict[int, Chem.Mol]
MonomerGraphMapping = ty.Dict[int, ty.Tuple[int, str]]
Reaction = ChemicalReaction 

class MolecularPattern:
    def __init__(self, name: str, core: bool, smarts: str) -> None:
        """
        Create a molecular pattern.

        Parameters
        ----------
        name : str
            Name of molecular pattern.
        core : bool
            Whether or not molecular pattern is a core unit.
        smarts : str
            SMARTS string of molecular pattern.
        
        Returns
        -------
        None
        """
        self.name = name 
        self.core = core
        self.smarts = smarts 
        self.compiled = self._compile_pattern(self.smarts)

    def is_core(self) -> bool:
        """
        Check if molecular pattern is a core unit.
        
        Parameters
        ----------
        None

        Returns
        -------
        is_core : bool
            Whether or not molecular pattern is a core unit.
        """
        return self.core

    def _compile_pattern(self, smarts: str) -> Chem.Mol:
        """
        Compile a molecular pattern.
        
        Parameters
        ----------
        smarts : str
            SMARTS string of molecular pattern.
        
        Returns
        -------
        compiled : Chem.Mol
            Compiled molecular pattern.
        """
        return Chem.MolFromSmarts(smarts)
    
class ReactionRule: 
    def __init__(self, name: str, smirks: str) -> None:
        """
        Create a reaction rule.
        
        Parameters
        ----------
        name : str
            Name of reaction rule.
        smirks : str
            SMIRKS string of reaction rule.
        
        Returns
        -------
        None
        """
        self.name = name 
        self.backward_smirks = smirks 
        self.forward_smirks = self._reverse_smirks(self.backward_smirks)
        self.backward_compiled = self._compile_rule(self.backward_smirks)
        self.forward_compiled = self._compile_rule(self.forward_smirks)

    def _reverse_smirks(self, smirks: str) -> str:
        """
        Reverse a SMIRKS string.
        
        Parameters
        ----------
        smirks : str
            SMIRKS string to reverse.
        
        Returns
        -------
        reversed_smirks : str
            Reversed SMIRKS string.
        """
        return ">>".join(reversed(smirks.split(">>")))

    def _compile_rule(self, smirks: str) -> Reaction:
        """
        Compile a reaction rule.
        
        Parameters
        ----------
        smirks : str
            SMIRKS string of reaction rule.
        
        Returns
        -------
        compiled : Reaction
            Compiled reaction rule.
        """
        return ReactionFromSmarts(smirks)
    
class Molecule:
    def __init__(self, name: str, smiles: str) -> None:
        """
        Create a molecule.
        
        Parameters
        ----------
        name : str
            Name of molecule.
        smiles : str
            SMILES string of molecule.
        
        Returns
        -------
        None
        """
        self.name = name 
        self.smiles = smiles 
        self.compiled = self._compile_molecule(self.smiles)

    def _compile_molecule(self, smiles: str) -> Chem.Mol:
        """
        Compile a molecule.
        
        Parameters
        ----------
        smiles : str
            SMILES string of molecule.
        
        Returns
        -------
        compiled : Chem.Mol
            Compiled molecule.
        """
        return Chem.MolFromSmiles(smiles)
    
    def apply_rules(
        self,
        reactions: ty.List[ReactionRule]
    ) -> ty.Tuple[Chem.Mol, Tree, ReactionTreeMapping]:
        """
        Apply reaction rules to a molecule.

        Parameters
        ----------
        reactions : ty.List[ReactionRule]
            List of reaction rules.
        
        Returns
        -------
        substrate : Chem.Mol
            Substrate.
        tree : Tree
            Reaction tree.
        mapping : ReactionTreeMapping
            Reaction tree mapping.
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

    Parameters
    ----------
    mol : Chem.Mol
        Molecule.
    radius : int
        Fingerprint radius.
    num_bits : int
        Number of bits in fingerprint.
    
    Returns
    -------
    fp_arr : np.array
        Fingerprint.    
    """
    fp_arr = np.zeros((0,), dtype=np.int8)
    fp_vec = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=num_bits)

    DataStructs.ConvertToNumpyArray(fp_vec, fp_arr)

    return fp_arr

def tanimoto_similarity(fp1: np.array, fp2: np.array) -> float:
    """
    Calculate Tanimoto similarity between two fingerprints.

    Parameters
    ----------
    fp1 : np.array
        First fingerprint.
    fp2 : np.array
        Second fingerprint.
    
    Returns
    -------
    similarity : float
        Tanimoto similarity.
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

    Parameters
    ----------
    mol : Chem.Mol
        Molecule.
    N : int
        Number of atoms in original parent substrate.
    radius : int
        Fingerprint radius.
    num_bits : int
        Number of bits in fingerprint.

    Returns
    -------
    hash : int
        Hash of reaction product.
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
    
    Parameters
    ----------
    mol : Chem.Mol
        Molecule.
    monomers : ty.List[MolecularPattern]
        List of monomers.
    
    Returns
    -------
    name : str
        Name of monomer.
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

    Parameters
    ----------
    mol : Chem.Mol
        Molecule.
    identified : ty.List[ty.Tuple[int, str]]
        List of identified monomers.
    mapping :  ReactionTreeMapping
        Mapping of reaction products to molecules.
    
    Returns
    -------
    selected_subsets : ty.List[ty.Tuple[int, str]]
        List of selected subsets, where each subset is a tuple of the form
        (node, node_id).
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