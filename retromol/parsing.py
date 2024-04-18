"""This module contains functions for parsing RetroMol Result objects into 
JSON and vice versa.
"""
import json
import typing as ty
from dataclasses import dataclass

from rdkit import Chem

from retromol.chem import Molecule, MolecularPattern, ReactionRule
from retromol.graph import reaction_tree_to_digraph, reaction_tree_to_monomer_graph

@dataclass
class Result:
    """Dataclass for storing the results of a RetroMol run.
    
    :param identifier: The identifier of the molecule.
    :type identifier: str
    :param mol: The molecule.
    :type mol: Chem.Mol
    :param success: Whether the RetroMol run was successful.
    :type success: bool
    :param score: The score of the molecule.
    :type score: int
    :param reaction_tree: The reaction tree.
    :type reaction_tree: ty.Dict[int, ty.List[int]]
    :param applied_reactions: The applied reactions.
    :type applied_reactions: ty.List[str]
    :param monomer_graph: The monomer graph.
    :type monomer_graph: ty.Dict[str, ty.Any]
    """
    identifier: str
    mol: Chem.Mol
    success: bool
    score: int = None
    reaction_tree: ty.Dict[int, ty.List[int]] = None
    applied_reactions: ty.List[str] = None
    monomer_graph: ty.Dict[str, ty.Any] = None

    def to_json(self, indent: int = 4) -> str:
        """Convert the result to JSON.
        
        :param indent: The number of spaces to indent the JSON.
        :type indent: int
        :return: The result as JSON.
        :rtype: str
        """
        return json.dumps({
            "identifier": self.identifier,
            "smiles": Chem.MolToSmiles(self.mol),
            "success": "true" if self.success else "false",
            "score": self.score,
            "reaction_tree": self.reaction_tree,
            "applied_reactions": self.applied_reactions,
            "monomer_graph": self.monomer_graph
        }, indent=indent)

    @classmethod
    def from_json(cls, path: str) -> "Result":
        """Create a Result object from a JSON file.
        
        :param path: The path to the JSON file.
        :type path: str
        :return: The Result object.
        :rtype: Result
        """
        with open(path, "r", encoding="utf-8") as fo:
            data = json.load(fo)

        reaction_tree = data["reaction_tree"]
        if reaction_tree is not None:
            reaction_tree = {int(k): v for k, v in reaction_tree.items()}

        monomer_graph = data["monomer_graph"]
        if monomer_graph is not None:
            monomer_graph = {int(k): v for k, v in monomer_graph.items()}

        return Result(
            identifier=data["identifier"],
            mol=Chem.MolFromSmiles(data["smiles"]),
            success=True if data["success"] == "true" else False,
            score=data["score"],
            reaction_tree=reaction_tree,
            applied_reactions=data["applied_reactions"],
            monomer_graph=monomer_graph,
        )

    def has_identified_monomers(self) -> bool:
        """Check if any monomers have been identified.
        
        :return: Whether any monomers have been identified.
        :rtype: bool
        """
        return any([x["identity"] is not None for _, x in self.monomer_graph.items()])

def parse_reaction_rules(src: str) -> ty.List[ReactionRule]:
    """Parse reaction rules from JSON.
    
    :param src: The JSON source.
    :type src: str
    :return: The reaction rules.
    :rtype: ty.List[ReactionRule]
    """
    data = json.loads(src)

    reaction_rules = []
    for item in data:
        try:
            reaction_rule = ReactionRule(item["identifier"], item["reaction_patterns"], item["properties"])
        except Exception as err:
            msg = f"{err}\nError parsing reaction rule:\n{item}"
            raise Exception(msg)

        reaction_rules.append(reaction_rule)

    return reaction_rules

def parse_molecular_patterns(src: str) -> ty.List[MolecularPattern]:
    """Parse molecular patterns from JSON.
    
    :param src: The JSON source.
    :type src: str
    :return: The molecular patterns.
    :rtype: ty.List[MolecularPattern]
    """
    data = json.loads(src)

    molecular_patterns = []
    for item in data:
        try:
            molecular_pattern = MolecularPattern(item["type"], item["patterns"], item["properties"])
        except Exception as err:
            msg = f"{err}\nError parsing molecular pattern:\n{item}"
            raise Exception(msg)
        molecular_patterns.append(molecular_pattern)

    return molecular_patterns

def parse_mol(
    mol: Molecule,
    reactions: ty.List[ReactionRule],
    monomers: ty.List[MolecularPattern],
) -> Result:
    """Parse a molecule with RetroMol.
    
    :param mol: The molecule.
    :type mol: Molecule
    :param reactions: The reaction rules.
    :type reactions: ty.List[ReactionRule]
    :param monomers: The molecular patterns.
    :type monomers: ty.List[MolecularPattern]
    :return: The result.
    :rtype: Result
    """
    name = mol.name
    reactant, reaction_tree, reaction_mapping = mol.apply_rules(reactions)
    applied_reactions = list(set([reaction for _, reactions in reaction_tree.items() for reaction in reactions]))
    reaction_tree = reaction_tree_to_digraph(reaction_tree)
    monomer_graph, monomer_mapping = reaction_tree_to_monomer_graph(mol, reaction_tree, reaction_mapping, monomers)

    # Reformat reaction_tree and reaction_mapping into one graph.
    new_reaction_tree = {}
    for parent in reaction_tree.nodes:
        children = list(reaction_tree.successors(parent))

        mol = reaction_mapping[parent]
        for atom in mol.GetAtoms():
            amn = atom.GetIsotope()
            atom.SetIsotope(0)
            atom.SetAtomMapNum(amn)
        smiles = Chem.MolToSmiles(mol)

        new_reaction_tree[parent] = {
            "smiles": smiles, 
            "children": children
        }

    # Reformat monomer_graph and monomer_mapping into one graph.
    new_monomer_graph = {}

    for reaction_tree_node, (monomer_graph_node, identity) in monomer_mapping.items():
        new_monomer_graph[monomer_graph_node] = {
            "reaction_tree_id": reaction_tree_node,
            "identity": identity,
            "neighbors": []
        }

    for node in monomer_graph:
        neighbors = list(monomer_graph[node])

        if node in new_monomer_graph:
            new_monomer_graph[node]["neighbors"] = neighbors
        else:
            new_monomer_graph[node] = {
                "reaction_tree_id": None,
                "identity": None,
                "neighbors": neighbors
            }

    # Create result object.
    return Result(
        identifier=name,
        mol=reactant,
        success=True,
        score=None,
        reaction_tree=new_reaction_tree,
        applied_reactions=applied_reactions,
        monomer_graph=new_monomer_graph,
    )
