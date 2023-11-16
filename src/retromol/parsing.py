import json 
import typing as ty 
from collections import Counter 
from dataclasses import dataclass

import matplotlib.pyplot as plt
import networkx as nx
from rdkit import Chem 

from retromol.chem import (
    Molecule, 
    MolecularPattern,
    ReactionRule,
    ReactionTreeMapping,
    MonomerGraphMapping,
    get_2d_coordinatates,
    mol_to_fingerprint,
    tanimoto_similarity
)
from retromol.graph import (
    reaction_tree_to_digraph,
    reaction_tree_to_monomer_graph,
    resolve_biosynthetic_sequence
)

@dataclass
class BiosyntheticRecipe:
    graph: nx.DiGraph
    coords: ty.Dict[int, ty.Tuple[float, float]]
    labels: ty.Dict[int, str]
    edge_labels: ty.Dict[ty.Tuple[int, int], str]
    node_color: ty.List[str]
    node_size: ty.List[int]

    def visualize(self, path: ty.Optional[str]) -> None:
        """
        Visualize biosynthetic recipe.

        Parameters
        ----------
        path : ty.Optional[str]
            Path to output png file.
        
        Returns
        -------
        None
        """
        nx.draw(self.graph, self.coords, node_color=self.node_color, node_size=self.node_size)
        nx.draw_networkx_labels(self.graph, self.coords, labels=self.labels)
        nx.draw_networkx_edge_labels(self.graph, self.coords, edge_labels=self.edge_labels, rotate=False)

        # Create legend.
        legend_elements = [
            plt.Line2D([0], [0], marker="o", color="w", label="Root", markerfacecolor="green", markersize=15),
            plt.Line2D([0], [0], marker="o", color="w", label="Core biosynthetic unit", markerfacecolor="red", markersize=15),
        ]
        plt.legend(handles=legend_elements)

        if path is None:
            plt.show()
        else:
            plt.savefig(path, dpi=300)
            plt.close()

def generate_structures(
    recipe: BiosyntheticRecipe,
    reactions: ty.List[ReactionRule],
    motifs: ty.List[MolecularPattern]
) -> ty.List[Chem.Mol]:
    """
    Generate structures.

    Parameters
    ----------
    recipe : BiosyntheticRecipe
        Biosynthetic recipe.
    reactions : ty.List[ReactionRule]
        List of reaction rules.
    motifs : ty.List[MolecularPattern]
        List of motifs.
        
    Returns 
    -------
    mols : ty.List[Chem.Mol]
        List of molecules.
    """
    print("Not yet implemented")  
    return []

def get_edge_labels(
    result: "Result", 
    root: int, 
    source_id: int, 
    target_id: int
) -> ty.List[str]:
    """
    Get edge label.
    
    Parameters
    ----------
    result : Result
        Result object.
    root : int
        Root node ID.
    source_id : int
        Source node ID.
    target_id : int
        Target node ID.
    
    Returns
    -------
    edge_labels : ty.List[str]
        Edge labels.
    """
    monomer_to_reaction_tree_id = {i[0]: j for j, i in result.monomer_mapping.items()}
    reaction_tree_source_id = monomer_to_reaction_tree_id[source_id]
    if target_id != root: # For start of sequence.
        reaction_tree_target_id = monomer_to_reaction_tree_id[target_id]
    else:
        reaction_tree_target_id = root

    # Get closest common ancestor in between source and target in reaction tree.
    source_path = nx.shortest_path(result.reaction_tree, root, reaction_tree_source_id)
    target_path = nx.shortest_path(result.reaction_tree, root, reaction_tree_target_id)
    common_ancestor = None
    for node in source_path:
        if node in target_path:
            common_ancestor = node
            break
    
    # Get reaction labels for all edges between common ancestor and source.
    source_path = nx.shortest_path(result.reaction_tree, common_ancestor, reaction_tree_source_id)
    source_path_edges = list(zip(source_path[:-1], source_path[1:]))
    source_path_edge_labels = [
        result.reaction_tree.edges[edge]["reaction"]
        for edge in source_path_edges
    ]
    
    # Get reaction labels for all edges between common ancestor and target.
    target_path = nx.shortest_path(result.reaction_tree, common_ancestor, reaction_tree_target_id)
    target_path_edges = list(zip(target_path[:-1], target_path[1:]))
    target_path_edge_labels = [
        result.reaction_tree.edges[edge]["reaction"]
        for edge in target_path_edges
    ]

    # Combine sets.
    # edge_labels = list(source_path_edge_labels | target_path_edge_labels)
    edge_labels = source_path_edge_labels

    # Sort labels.
    edge_labels = sorted(edge_labels)
    
    return list(edge_labels)

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
    
    def get_biosynthetic_recipe(self) -> BiosyntheticRecipe:
        """
        Get biosynthetic recipe.

        Parameters
        ----------
        None
        
        Returns
        -------
        graph : nx.DiGraph
            Graph.
        """
        # Only print first one. 
        seqs = self.biosynthetic_seq
        if len(seqs) != 1:
            raise ValueError(f"Expected one biosynthetic sequence, got {len(seqs)}.")
        seq = seqs[0]

        # Get root of reaction graph.
        roots = [
            node 
            for node in self.reaction_tree.nodes 
            if self.reaction_tree.in_degree(node) == 0
        ]
        if len(roots) != 1:
            raise ValueError(f"Expected one root node, got {len(roots)}.")
        root = roots[0]

        # Add root to seq.
        seq = [(root, "root")] + seq

        edge_labels = {}
        for i in range(len(seq) - 1):
            rev_seq = seq[::-1]
            source_id, source_name = rev_seq[i]
            target_id, target_name = rev_seq[i + 1]
            labels = get_edge_labels(self, root, source_id, target_id)
            label_counter = Counter(labels)
            label_counter_str = " / ".join([f"{k} ({v})" for k, v in label_counter.items()])
            edge_labels[(source_id, target_id)] = label_counter_str

        graph = nx.DiGraph()
        coords = {}
        labels = {}
        node_color = []
        node_size = []
        for i in range(len(seq) - 1):
            source_id, source_name = seq[i]
            target_id, target_name = seq[i + 1]
            graph.add_node(source_id)
            graph.add_node(target_id)
            graph.add_edge(source_id, target_id)
            
            coords[source_id] = (0, -i)
            if source_name == "root":
                node_color.append("green")
                node_size.append(100)
            else:
                node_color.append("red")
                node_size.append(200)
                labels[source_id] = source_name

            if i == len(seq) - 2:
                coords[target_id] = (0, -i - 1)
                labels[target_id] = target_name
                node_color.append("red")
                node_size.append(200)

        return BiosyntheticRecipe(
            graph,
            coords,
            labels,
            edge_labels,
            node_color,
            node_size
        )
    
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

        self.name = data["name"]
        self.substrate = Chem.MolFromSmiles(data["smiles"])
        self.score = data["score"]

        monomer_graph_json = data["monomer_graph"]
        monomer_graph = nx.readwrite.json_graph.node_link_graph(monomer_graph_json)

        monomer_mapping = data["monomer_mapping"]
        monomer_mapping = {int(k): (int(v), s) for k, (v, s) in monomer_mapping.items()}

        reaction_tree_json = data["reaction_tree"]
        reaction_tree = nx.readwrite.json_graph.node_link_graph(reaction_tree_json)

        reaction_mapping = data["reaction_mapping"]
        reaction_mapping = {int(k): Chem.MolFromSmiles(m) for k, m in reaction_mapping.items()}

        biosynthetic_seq = data["biosynthetic_seq"]

        result = Result(
            self.name,
            self.substrate,
            self.score,
            reaction_tree,
            reaction_mapping,
            monomer_graph,
            monomer_mapping,
            biosynthetic_seq
        )

        return result
    
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

def parse_molecular_patterns(path: str, as_smiles: bool = False) -> ty.List[MolecularPattern]:   
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
            molecular_pattern = MolecularPattern(item["name"], item["smarts"], as_smiles=as_smiles)
        except Exception as err:
            msg = f"{err}\nError parsing molecular pattern:\n{item}"
            raise Exception(msg)
        
        molecular_patterns.append(molecular_pattern)

    return molecular_patterns

def parse_mol(
    mol: Molecule, 
    reactions: ty.List[ReactionRule],
    core_units: ty.List[MolecularPattern],
    other_units: ty.List[MolecularPattern],
) -> Result:
    """
    Parse molecule.

    Parameters
    ----------
    mol : Molecule
        Molecule.
    reactions : ty.List[ReactionRule]
        List of reaction rules.
    core_units : ty.List[MolecularPattern]
        List of motif units.
    other_units : ty.List[MolecularPattern]
        List of other monomer units.
    
    Returns
    -------
    result : Result
        Result object.
    """
    substrate, reaction_tree, reaction_mapping = mol.apply_rules(reactions)
    reaction_tree = reaction_tree_to_digraph(reaction_tree)

    # # Print all reactions.
    # print(len(reaction_tree))
    # reactions = set()
    # for edge in reaction_tree.edges:
    #     reaction = reaction_tree.edges[edge]["reaction"]
    #     reactions.add(reaction)
    # print(reactions)

    # exit()

    monomers = core_units + other_units
    monomer_graph, monomer_mapping = reaction_tree_to_monomer_graph(mol, reaction_tree, reaction_mapping, monomers)
    biosynthetic_seq = resolve_biosynthetic_sequence(reaction_tree, reaction_mapping, monomer_graph, monomer_mapping, core_units)
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

def visualize_monomer_graph(data: Result, path: ty.Optional[str]) -> None:
    """
    Visualize monomer graph.
    
    Parameters
    ----------
    data : Result
        Result object.
    path : ty.Optional[str]
        Path to output png file.
    
    Returns
    -------
    None
    """
    # Get 2D coordinates of atoms in molecule. The atom tracking numbers are 
    # also used to identify nodes in the monomer graph.
    pos = get_2d_coordinatates(data.substrate)

    # Color nodes based on if node resembles monomer.
    monomer_ids = [v[0] for v in data.monomer_mapping.values()]
    node_color = []
    node_size = []
    for node in data.monomer_graph.nodes:
        if node in monomer_ids:
            node_color.append("red")
            node_size.append(200)
        else:
            node_color.append("blue")
            node_size.append(50)
    nx.draw(data.monomer_graph, pos, node_color=node_color, node_size=node_size)

    # Draw identity label when identity of node is known.
    labels = {}
    for _, (node, identity) in data.monomer_mapping.items():
        labels[node] = identity
    nx.draw_networkx_labels(data.monomer_graph, pos, labels=labels)

    # Create legend.
    legend_elements = [
        plt.Line2D([0], [0], marker="o", color="w", label="Atom", markerfacecolor="blue", markersize=15),
        plt.Line2D([0], [0], marker="o", color="w", label="Monomer", markerfacecolor="red", markersize=15),
    ]
    plt.legend(handles=legend_elements)

    if path is None:
        plt.show()
    else:
        plt.savefig(path, dpi=300)
        plt.close()