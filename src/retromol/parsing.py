import json 
import typing as ty 
from dataclasses import dataclass

from rdkit import Chem 

from retromol.chem import Molecule, MolecularPattern, ReactionRule
from retromol.graph import reaction_tree_to_digraph, reaction_tree_to_monomer_graph

@dataclass 
class Result:
    identifier: str
    substrate: Chem.Mol
    success: bool
    score: int = None
    reaction_tree: ty.Dict[int, ty.List[int]] = None 
    monomer_graph: ty.Dict[str, ty.Any] = None

    def to_json(self, indent: int = 4) -> str:
        return json.dumps({
            "identifier": self.identifier,
            "smiles": Chem.MolToSmiles(self.substrate),
            "success": "true" if self.success else "false",
            "score": self.score,
            "reaction_tree": self.reaction_tree,
            "monomer_graph": self.monomer_graph
        }, indent=indent)
    
    @classmethod
    def from_json(self, path: str) -> "Result":
        with open(path, "r") as file_open:
            data = json.load(file_open)

        return Result(
            identifier=data["identifier"],
            substrate=Chem.MolFromSmiles(data["smiles"]),
            succes=True if data["success"] == "true" else False,
            score=data["score"],
            reaction_tree=data["reaction_tree"],
            monomer_graph=data["monomer_graph"],
        )
    
def parse_reaction_rules(src: str) -> ty.List[ReactionRule]:
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
    name = mol.name
    substrate, reaction_tree, reaction_mapping = mol.apply_rules(reactions)
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
        substrate=substrate, 
        success=True,
        score=None,
        reaction_tree=new_reaction_tree,
        monomer_graph=new_monomer_graph,
    )