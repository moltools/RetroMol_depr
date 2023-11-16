#!/usr/bin/env python3
import argparse 
import os

from retromol.chem import draw_molecule
from retromol.parsing import Result, visualize_monomer_graph

def cli() -> argparse.Namespace:   
    """
    Command-line interface.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="Path to JSON file containing monomer graph.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to output dir.")
    return parser.parse_args()

def main() -> None:
    """
    Driver code.
    """
    args = cli()
    result = Result.from_json(args.input)

    out_molecule_img = os.path.join(args.output, "molecule.png")
    out_monomer_graph = os.path.join(args.output, "monomer_graph.png")
    out_biosynthetic_graph = os.path.join(args.output, "biosynthetic_graph.png")

    draw_molecule(result.substrate, out_molecule_img)
    visualize_monomer_graph(result, out_monomer_graph) 
    recipe = result.get_biosynthetic_recipe()
    recipe.visualize(out_biosynthetic_graph)

    exit(0)

if __name__ == "__main__":
    main()
