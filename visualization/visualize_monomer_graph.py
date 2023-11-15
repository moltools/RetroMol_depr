#!/usr/bin/env python3
import argparse 
import json
import typing as ty 

import networkx as nx 
import matplotlib.pyplot as plt 
from rdkit import Chem 
from rdkit.Chem import AllChem 

def cli() -> argparse.Namespace:   
    """
    Command-line interface.
    
    Parameters
    ----------
    None

    Returns
    ------- 
    args : argparse.Namespace
        Namespace object containing the parsed arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Path to JSON file containing monomer graph."
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Path to output png file."
    )
    return parser.parse_args()

def read_monomer_graph(path: str) -> ty.Tuple[Chem.Mol, nx.Graph, ty.Dict[int, int]]:
    """
    Read monomer graph from JSON file.
    
    Parameters
    ----------
    path : str
        Path to JSON file containing monomer graph.
    
    Returns
    -------
    mol : Chem.Mol
        Substrate.
    monomer_graph : nx.Graph
        Monomer graph.
    monomer_mapping : ty.Dict[int, int]
        Mapping from reaction graph node to monomer graph node.
    """
    with open(path, "r") as file_open:
        data = json.load(file_open)

    mol = Chem.MolFromSmiles(data["smiles"])

    monomer_graph_json = data["monomer_graph"]
    monomer_graph = nx.readwrite.json_graph.node_link_graph(monomer_graph_json)

    monomer_mapping = data["monomer_mapping"]
    monomer_mapping = {int(k): (int(v), s) for k, (v, s) in monomer_mapping.items()}
    
    return mol, monomer_graph, monomer_mapping

def get_2d_coordinatates(mol: Chem.Mol) -> ty.Dict[int, ty.Tuple[float, float]]:
    """
    Get 2D coordinates of atoms in molecule.
    
    Parameters
    ----------
    mol : Chem.Mol
        Molecule.
    
    Returns
    -------
    coordinates : ty.Dict[int, ty.Tuple[float, float]]
        Mapping from atom index to 2D coordinate.
    """
    AllChem.Compute2DCoords(mol)

    coordinates = {}
    for atom in mol.GetAtoms():
        atom_idx = atom.GetIdx()
        position = mol.GetConformer(0).GetAtomPosition(atom_idx)
        amn = atom.GetIsotope() # Stored atom tracking num as isotope prop.
        coordinates[amn] = (position.x, position.y)

    return coordinates

def visualize_monomer_graph(
    mol: Chem.Mol,
    monomer_graph: nx.Graph, 
    monomer_mapping: ty.Dict[int, int],
    path: str
) -> None:
    """
    Visualize monomer graph.
    
    Parameters
    ----------
    mol : Chem.Mol
        Substrate.
    monomer_graph : nx.Graph
        Monomer graph.
    monomer_mapping : ty.Dict[int, int]
        Mapping from reaction graph node to monomer graph node.
    path : str
        Path to output png file.
    
    Returns
    -------
    None
    """
    # Get 2D coordinates of atoms in molecule. The atom tracking numbers are 
    # also used to identify nodes in the monomer graph.
    pos = get_2d_coordinatates(mol)

    # Color nodes based on if node resembles monomer.
    monomer_ids = [v[0] for v in monomer_mapping.values()]
    node_color = []
    node_size = []
    for node in monomer_graph.nodes:
        if node in monomer_ids:
            node_color.append("red")
            node_size.append(200)
        else:
            node_color.append("blue")
            node_size.append(50)
    nx.draw(monomer_graph, pos, node_color=node_color, node_size=node_size)

    # Draw identity label when identity of node is known.
    labels = {}
    for _, (node, identity) in monomer_mapping.items():
        labels[node] = identity
    nx.draw_networkx_labels(monomer_graph, pos, labels=labels)

    # Create legend.
    legend_elements = [
        plt.Line2D([0], [0], marker="o", color="w", label="Atom", markerfacecolor="blue", markersize=15),
        plt.Line2D([0], [0], marker="o", color="w", label="Monomer", markerfacecolor="red", markersize=15),
    ]
    plt.legend(handles=legend_elements)
    plt.savefig(path, dpi=300)
    plt.close()

def main() -> None:
    """
    Driver code.
    """
    args = cli()
    mol, graph, mapping = read_monomer_graph(args.input)
    visualize_monomer_graph(mol, graph, mapping, args.output) 
    exit(0)

if __name__ == "__main__":
    main()
