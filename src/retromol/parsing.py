import json 
import typing as ty 
from dataclasses import dataclass

import networkx as nx
from rdkit import Chem 

from retromol.chem import Molecule, MolecularPattern, ReactionRule, ReactionTreeMapping, MonomerGraphMapping
from retromol.graph import reaction_tree_to_digraph, reaction_tree_to_monomer_graph, resolve_biosynthetic_sequence

@dataclass 
class Result:
    name: str
    substrate: Chem.Mol
    success: bool
    score: int = None
    reaction_tree: nx.DiGraph = None
    reaction_mapping: ReactionTreeMapping = None
    monomer_graph: nx.Graph = None
    monomer_mapping: MonomerGraphMapping = None
    biosynthetic_seq: ty.List[ty.Tuple[int, str]] = None

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
        if self.success:
            reaction_tree_json = nx.readwrite.json_graph.node_link_data(self.reaction_tree)
            reaction_mapping = self._reaction_mapping_to_json(self.reaction_mapping)
            monomer_graph_json = nx.readwrite.json_graph.node_link_data(self.monomer_graph)
        
            # Reformat links so that they are always JSON serializable.
            monomer_graph_links = []
            for link in monomer_graph_json["links"]:
                new_link = {"source": link["source"], "target": link["target"]}
                monomer_graph_links.append(new_link)
            monomer_graph_json["links"] = monomer_graph_links

        else:
            reaction_tree_json = None
            reaction_mapping = None
            monomer_graph_json = None

        data = {
            "name": self.name,
            "smiles": Chem.MolToSmiles(self.substrate),
            "success": "true" if self.success else "false",
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
    
    @classmethod
    def from_json(self, path: str) -> "Result":
        """
        Convert JSON string to result.

        Parameters
        ----------
        path : str
            Path to JSON file.
        
        Returns
        -------
        result : Result
            Result object.
        """
        with open(path, "r") as file_open:
            data = json.load(file_open)

        name = data["name"]
        substrate = Chem.MolFromSmiles(data["smiles"])
        success = True if data["success"] == "true" else False

        if success:
            score = data["score"]

            monomer_graph_json = data["monomer_graph"]
            monomer_graph = nx.readwrite.json_graph.node_link_graph(monomer_graph_json)

            monomer_mapping = data["monomer_mapping"]
            monomer_mapping = {int(k): (int(v), s) for k, (v, s) in monomer_mapping.items()}

            reaction_tree_json = data["reaction_tree"]
            reaction_tree = nx.readwrite.json_graph.node_link_graph(reaction_tree_json)

            reaction_mapping = data["reaction_mapping"]
            reaction_mapping = {int(k): Chem.MolFromSmiles(m) for k, m in reaction_mapping.items()}

            biosynthetic_seq = data["biosynthetic_seq"]
        
        else:
            score = None
            reaction_tree = None
            reaction_mapping = None
            monomer_graph = None
            monomer_mapping = None
            biosynthetic_seq = None

        result = Result(
            name,
            substrate,
            success,
            score,
            reaction_tree,
            reaction_mapping,
            monomer_graph,
            monomer_mapping,
            biosynthetic_seq
        )

        return result
    
def parse_reaction_rules(src: str) -> ty.List[ReactionRule]:
    """
    Parse reaction rules from file.
    
    Parameters
    ----------
    src : str 
        String in JSON format describing reaction rules.
    
    Returns
    -------
    reaction_rules : ty.List[ReactionRule]
        List of reaction rules.
    """
    reaction_rules = []
    
    data = json.loads(src)
    
    for item in data:
        try:
            reaction_rule = ReactionRule(item["name"], item["smirks"])
        except Exception as err:
            msg = f"{err}\nError parsing reaction rule:\n{item}"
            raise Exception(msg)

        reaction_rules.append(reaction_rule)

    return reaction_rules

def parse_molecular_patterns(src: str) -> ty.List[MolecularPattern]:   
    """
    Parse molecular patterns from file.
    
    Parameters
    ----------
    src : str
        String in JSON format describing molecular patterns.
    
    Returns
    -------
    molecular_patterns : ty.List[MolecularPattern]
        List of molecular patterns.
    """
    molecular_patterns = []

    data = json.loads(src)
    
    for item in data:
        try:
            molecular_pattern = MolecularPattern(item["name"], item["core"], item["smarts"])
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
    """
    Parse molecule.

    Parameters
    ----------
    mol : Molecule
        Molecule.
    reactions : ty.List[ReactionRule]
        List of reaction rules.
    monomers : ty.List[MolecularPattern]
        List of monomer units.
    
    Returns
    -------
    result : Result
        Result object.
    """
    substrate, reaction_tree, reaction_mapping = mol.apply_rules(reactions)
    reaction_tree = reaction_tree_to_digraph(reaction_tree)

    core_monomers = [m for m in monomers if m.is_core()]
    monomer_graph, monomer_mapping = reaction_tree_to_monomer_graph(mol, reaction_tree, reaction_mapping, monomers)
    biosynthetic_seq = resolve_biosynthetic_sequence(reaction_tree, reaction_mapping, monomer_graph, monomer_mapping, core_monomers)
    score = len(monomer_graph.nodes) - len(monomer_mapping)
    
    result = Result(
        mol.name,
        substrate, 
        True,
        score,
        reaction_tree,
        reaction_mapping,
        monomer_graph,
        monomer_mapping,
        biosynthetic_seq
    )

    return result