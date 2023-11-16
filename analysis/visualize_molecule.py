#!/usr/bin/env python3
import argparse

from rdkit import Chem

from retromol.chem import draw_molecule

def cli() -> argparse.Namespace:
    """
    Command-line interface.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="Input SMILES.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to output path.")
    return parser.parse_args()

def main() -> None:
    """
    Driver code.
    """
    args = cli()
    draw_molecule(Chem.MolFromSmiles(args.input), args.output)
    exit(0)

if __name__ == "__main__":
    main()