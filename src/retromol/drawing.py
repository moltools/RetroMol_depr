import typing as ty 
from copy import deepcopy

import matplotlib.pyplot as plt
import networkx as nx
from rdkit import Chem 
from rdkit.Chem import AllChem, Draw

from retromol.parsing import Result

def get_2d_coordinatates(mol: Chem.Mol) -> ty.Dict[int, ty.Tuple[float, float]]:
    """
    NOTE: set desired atom mapping as atom isotope number.
    """
    AllChem.Compute2DCoords(mol)

    coordinates = {}
    for atom in mol.GetAtoms():
        atom_idx = atom.GetIdx()
        position = mol.GetConformer(0).GetAtomPosition(atom_idx)
        amn = atom.GetIsotope() # Stored atom tracking num as isotope prop.
        coordinates[amn] = (position.x, position.y)

    return coordinates

def draw_molecule(mol: Chem.Mol, path: str) -> None:
    mol = deepcopy(mol) # Keep isotope numbers intact for atom mapping outside of function.

    for atom in mol.GetAtoms():
        atom.SetIsotope(0)

    AllChem.Compute2DCoords(mol)
    Draw.MolToImageFile(mol, path, size=(1000, 1000))

def visualize_monomer_graph(data: Result, path: ty.Optional[str]) -> None:
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