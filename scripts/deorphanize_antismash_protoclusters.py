#!/usr/bin/env python3
"""
Description:    This script finds nearest neighbors for protoclusters in the biosynthetic embedding space. 
Usage:          python3 find_nearest_neighbors.py -i <path/to/embedding> -m <path/to/embedder>
"""
import argparse
import joblib 

import numpy as np
from sklearn.neighbors import KDTree

from retromol_sequencing.fingerprint import get_biosynthetic_fingerprint

def cli() -> argparse.Namespace:
    """
    Command line interface.

    :returns: Command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="Path to input .npy embedding.")
    parser.add_argument("-n", "--labels", type=str, required=True, help="Path to input .txt embeddinglabels.")
    parser.add_argument("-m", "--embedder", type=str, required=True, help="Path to embedding model.")
    return parser.parse_args()

def main() -> None:
    """
    Main function.
    """
    args = cli()

    # Load embedding and embedder.
    embedding = np.load(args.input)
    labels = np.loadtxt(args.labels, dtype=str)
    embedder = joblib.load(args.embedder)

    # Deorphanization.
    seq = [""]
    query = embedder.transform([get_biosynthetic_fingerprint(seq)])
    
    # Find nearest neighbors with KDTree.
    tree = KDTree(embedding)
    dist, ind = tree.query(query, k=10)

    # Print results.
    print("Nearest neighbors:")
    for i, index in enumerate(ind[0]):
        print(f"\t{i+1}. {labels[index]} (distance: {round(dist[0][i], 3)})")

    exit(0)

if __name__ == "__main__":
    main()