#!/usr/bin/env python3
"""
Description:    Script to create Principal Moments of Inertia (PMI) plot of a set of molecular conformers.
                See https://pubs.acs.org/doi/full/10.1021/ci025599w# for more information on the PMI plot.
"""
import argparse 
import typing as ty

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from rdkit import Chem
from scipy.linalg import eig 

def cli() -> argparse.Namespace:
    """
    Command line interface.

    :returns: Command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="Input tsv file as 'name\tsmiles\n'.")
    parser.add_argument("--h", action="store_true", help="Whether the input file has a header.")
    parser.add_argument("-o", "--out", type=str, required=True, help="Output html file.")
    return parser.parse_args()

def calculate_center_of_mass(coords: np.ndarray, weights: ty.List[float]) -> np.ndarray:
    """
    Calculate center of mass of a set of coordinates.

    :param np.ndarray coords: Coordinates.
    :param ty.List[float] weights: Weights.
    :returns: Center of mass.
    :rtype: np.ndarray
    """
    return np.average(coords, axis=0, weights=weights)

def calculate_moments_of_inertia_and_axes(coords: np.ndarray, weights: ty.List[float]) -> np.ndarray:
    """
    Calculate principal moments of inertia and principal axes of a set of coordinates.

    :param np.ndarray coords: Coordinates.
    :param ty.List[float] weights: Weights.
    :returns: Principal moments of inertia and principal axes.
    :rtype: np.ndarray
    """
    # Check if weights is a list of length n.
    if not isinstance(weights, list):
        raise TypeError(f"weights must be a list, not {type(weights)}.")
    if len(weights) != coords.shape[0]:
        raise ValueError(f"weights must be a list of length {coords.shape[0]}, not {len(weights)}.")
    
    # Check if weights is a list of length n.
    if not isinstance(coords, np.ndarray):
        raise TypeError(f"coords must be a np.ndarray, not {type(coords)}.")
    if coords.shape[1] != 3:
        raise ValueError(f"coords must be a np.ndarray of shape (n, 3), not {coords.shape}.")
    
    # Calculate center of mass.
    center_of_mass = calculate_center_of_mass(coords, weights)
    
    # Translate center of mass.
    coords -= center_of_mass

    # Calculate moment of inertia tensor.
    tensor = np.zeros((3, 3))
    for i, (x, y, z) in enumerate(coords):
        tensor[0, 0] += weights[i] * (y**2 + z**2)
        tensor[1, 1] += weights[i] * (x**2 + z**2)
        tensor[2, 2] += weights[i] * (x**2 + y**2)
        tensor[0, 1] -= weights[i] * x * y
        tensor[0, 2] -= weights[i] * x * z
        tensor[1, 2] -= weights[i] * y * z
    tensor[1, 0] = tensor[0, 1]
    tensor[2, 0] = tensor[0, 2]
    tensor[2, 1] = tensor[1, 2]

    # Calcualte principal axes. 
    eigvals, eigvecs = eig(tensor)
    idx = eigvals.argsort()
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    return eigvals, eigvecs

def main() -> None:
    """
    Main function.
    """
    # Parse command line arguments.
    args = cli()

    # Parse mols from SDF file.
    mol_supplier = Chem.SDMolSupplier(args.input, removeHs=False)

    # Calculate principal moments of inertia.
    ratios = []
    for i, mol in enumerate(mol_supplier):
        try:
            # Try to get name.
            try:
                name = mol.GetProp("_Name")
            except:
                name = f"mol_{i}"

            # Get conformer.
            conf = mol.GetConformer()
            
            # Get positions and weights. 
            coords = conf.GetPositions()
            weights = [atom.GetMass() for atom in mol.GetAtoms()]

            # Calculate principal moments of inertia and principal axes.
            moments, _ = calculate_moments_of_inertia_and_axes(coords, weights)

            # Get normalized PMI ratios.
            i13 = moments[0] / moments[2]
            i23 = moments[1] / moments[2]

            # Append to list.
            ratios.append((name, float(i13), float(i23)))
        
        except Exception as e:
            print(f"Error: {e}")
            print(f"Skipping {name}.")
            continue

        print(f"{i}", end="\r")

    # Plot PMI ratios.
    labels = [ratio[0] for ratio in ratios]
    x = [ratio[1] for ratio in ratios]
    y = [ratio[2] for ratio in ratios]
    fig = px.scatter(x=x, y=y, hover_name=labels)
    fig.update_traces(marker=dict(symbol="triangle-up", size=15, color="blue"))

    triangle_x = [0.0, 0.5, 1.0, 0.0]
    triangle_y = [1.0, 0.5, 1.0, 1.0]
    fig.add_trace(go.Scatter(
        x=triangle_x, 
        y=triangle_y, 
        line=dict(color="black", width=2, dash="dash"),
        mode="lines",
        name="Shape envelope"    
    ))

    fig.add_annotation(text="Rod", x=-0.025, y=1.025, showarrow=False, font=dict(size=14, family="serif"))
    fig.add_annotation(text="Disc", x=0.5, y=0.475, showarrow=False, font=dict(size=14, family="serif"))
    fig.add_annotation(text="Sphere", x=1.025, y=1.025, showarrow=False, font=dict(size=14, family="serif"))

    fig.update_xaxes(
        title_text="$I_1/I_3$", 
        scaleanchor="y", 
        scaleratio=1, 
        tickfont=dict(size=16, family="serif"),
        linecolor="black",
    )
    fig.update_yaxes(
        title_text="$I_2/I_3$", 
        scaleanchor="x", 
        scaleratio=1, 
        tickfont=dict(size=16, family="serif"),
        linecolor="black",
    )
    fig.update_layout(
        font=dict(size=18, family="serif"), 
        title_text="Normalized principal moments of inertia",
        paper_bgcolor="white", 
        plot_bgcolor="white",
    )
    fig.write_html(
        args.out, 
        include_mathjax="cdn" # To save LaTeX annotations appropriately.
    )

    exit(0)

if __name__ == "__main__":
    main()