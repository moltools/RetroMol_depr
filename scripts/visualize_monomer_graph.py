#!/usr/bin/env python3
"""
Description:    Visualize monomer graph.
Usage:          python3 visualize_monomer_graph.py -i <path to RetroMol JSON result file> -o <path to output PNG>
"""
import argparse

from retromol.parsing import Result 
from retromol.drawing import visualize_monomer_graph

def cli() -> argparse.Namespace:
    """
    Command-line interface.

    :return:    Parsed command-line arguments.
    :rtype:     argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="Path to RetroMol JSON result file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to new PNG output path.")
    return parser.parse_args()

def main() -> None:
    """
    Driver code.
    """
    args = cli()
    result = Result.from_json(args.input)
    visualize_monomer_graph(result, args.output)
    exit(0)

if __name__ == "__main__":
    main()