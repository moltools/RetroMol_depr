import json 
import typing as ty 
from dataclasses import dataclass

import networkx as nx
from rdkit import Chem 

from retromol.chem import (
    Molecule, 
    MolecularPattern,
    ReactionRule,
    ReactionTreeMapping,
    MonomerGraphMapping
)
from retromol.graph import (
    reaction_tree_to_digraph,
    reaction_tree_to_monomer_graph,
    resolve_biosynthetic_sequence
)

@dataclass 
class Result:
    name: str
    substrate: Chem.Mol
    score: int
    reaction_tree: nx.DiGraph
    reaction_mapping: ReactionTreeMapping
    monomer_graph: nx.Graph
    monomer_mapping: MonomerGraphMapping
    biosynthetic_seq: ty.List[ty.Tuple[int, str]]

    def to_json(self, indent: int = 4) -> str:
        """
        Convert result to JSON string.

        Parameters
        ----------
        indent : int, optional
            Indentation level, by default 4.

        Returns
        -------
        data : str
            JSON string.
        """
        reaction_tree_json = nx.readwrite.json_graph.node_link_data(self.reaction_tree)
        reaction_mapping = self._reaction_mapping_to_json(self.reaction_mapping)
        monomer_graph_json = nx.readwrite.json_graph.node_link_data(self.monomer_graph)
        
        # Reformat links so that they are always JSON serializable.
        monomer_graph_links = []
        for link in monomer_graph_json["links"]:
            new_link = {"source": link["source"], "target": link["target"]}
            monomer_graph_links.append(new_link)
        monomer_graph_json["links"] = monomer_graph_links

        data = {
            "name": self.name,
            "smiles": Chem.MolToSmiles(self.substrate),
            "score": self.score,
            "reaction_tree": reaction_tree_json,
            "reaction_mapping": reaction_mapping,
            "monomer_graph": monomer_graph_json,
            "monomer_mapping": self.monomer_mapping,
            "biosynthetic_seq": self.biosynthetic_seq
        }
        return json.dumps(data, indent=indent)
    
    def _reaction_mapping_to_json(self, mapping: ReactionTreeMapping) -> str:
        """
        Convert mapping to JSON string.

        Parameters
        ----------
        mapping : Mapping
            Mapping.
        
        Returns
        -------
        data : str
            JSON string.
        """
        data = {}
    
        for node_id, node_mol in mapping.items():
            node_smiles = Chem.MolToSmiles(node_mol)
            data[node_id] = node_smiles

        return data
    
def parse_reaction_rules(path: str) -> ty.List[ReactionRule]:
    """
    Parse reaction rules from file.
    
    Parameters
    ----------
    path : str
        Path to file containing reaction rules.
    
    Returns
    -------
    reaction_rules : ty.List[ReactionRule]
        List of reaction rules.
    """
    reaction_rules = []
    
    with open(path, "r") as file_open:
        data = json.load(file_open)
    
    for item in data:
        try:
            reaction_rule = ReactionRule(item["name"], item["smirks"])
        except Exception as err:
            msg = f"{err}\nError parsing reaction rule:\n{item}"
            raise Exception(msg)

        reaction_rules.append(reaction_rule)

    return reaction_rules

def parse_molecular_patterns(path: str) -> ty.List[MolecularPattern]:   
    """
    Parse molecular patterns from file.
    
    Parameters
    ----------
    path : str
        Path to file containing molecular patterns.
    
    Returns
    -------
    molecular_patterns : ty.List[MolecularPattern]
        List of molecular patterns.
    """
    molecular_patterns = []

    with open(path, "r") as file_open:
        data = json.load(file_open)
    
    for item in data:
        try:
            molecular_pattern = MolecularPattern(item["name"], item["smarts"])
        except Exception as err:
            msg = f"{err}\nError parsing molecular pattern:\n{item}"
            raise Exception(msg)
        
        molecular_patterns.append(molecular_pattern)

    return molecular_patterns

def parse_mol(
    mol: Molecule, 
    reactions: ty.List[ReactionRule],
    motif_units: ty.List[MolecularPattern],
    starter_units: ty.List[MolecularPattern],
    tailoring_units: ty.List[MolecularPattern]
) -> Result:
    """
    Parse molecule.

    Parameters
    ----------
    mol : Molecule
        Molecule.
    reactions : ty.List[ReactionRule]
        List of reaction rules.
    motif_units : ty.List[MolecularPattern]
        List of motif units.
    starter_units : ty.List[MolecularPattern]
        List of starter units.
    tailoring_units : ty.List[MolecularPattern]
        List of tailoring units.
    
    Returns
    -------
    result : Result
        Result object.
    """
    substrate, reaction_tree, reaction_mapping = mol.apply_rules(reactions)
    reaction_tree = reaction_tree_to_digraph(reaction_tree)
    monomers = motif_units + starter_units + tailoring_units
    monomer_graph, monomer_mapping = reaction_tree_to_monomer_graph(mol, reaction_tree, reaction_mapping, monomers)
    biosynthetic_seq = resolve_biosynthetic_sequence(reaction_tree, reaction_mapping, monomer_graph, monomer_mapping, motif_units)
    score = len(monomer_graph.nodes) - len(monomer_mapping)
    
    result = Result(
        mol.name,
        substrate, 
        score,
        reaction_tree,
        reaction_mapping,
        monomer_graph,
        monomer_mapping,
        biosynthetic_seq
    )

    return result