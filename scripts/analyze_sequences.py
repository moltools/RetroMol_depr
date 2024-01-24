#!/usr/bin/env python3
"""
Description:    This script creates a primary sequence embedding space.
Usage:          python3 create_embedding.py -i <path/to/retromol/sequences/file> -o <path/to/output/dir>
"""
import argparse
import joblib 
import re
import typing as ty

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import umap 
from sklearn.cluster import HDBSCAN

from retromol_sequencing.fingerprint import get_biosynthetic_fingerprint, CompoundClassMapping
from retromol_sequencing.aligner import ModuleSequence, MultipleSequenceAlignment

def cli() -> argparse.Namespace:
    """
    Command line interface.

    :returns: Command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="Path to input .seq file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to output dir.")
    return parser.parse_args()

def plot_msa(path: str, msa: ty.List[ModuleSequence]) -> None:
    """
    Plot multiple sequence alignment.

    :param str path: Path to output file.
    :param ty.List[ModuleSequence] msa: Multiple sequence alignment.
    """
    # Get sequence labels.
    seq_labels = [seq.name for seq in msa]

    # Scale figure size with number and length of sequences.
    figsize_x = len(msa[0]._seq)
    figsize_y = len(msa)
    fig, ax = plt.subplots(figsize=(figsize_x, figsize_y))

    for i, seq in enumerate(msa):
        for j, item in enumerate(seq._seq):
            item = item[0]
            fontsize = 12
            if item == CompoundClassMapping.Undefined:
                facecolor, edgecolor = "white", "white"
                label = ""
            elif item == CompoundClassMapping.A:
                facecolor, edgecolor = "red", "black"
                label = item.name
            elif item == CompoundClassMapping.B:
                facecolor, edgecolor = "blue", "black"
                label = item.name
            elif item == CompoundClassMapping.C:
                facecolor, edgecolor = "green", "black"
                label = item.name
            elif item == CompoundClassMapping.D:
                facecolor, edgecolor = "orange", "black"
                label = item.name
            else:
                label = item.name
                facecolor, edgecolor = "purple", "black"
                fontsize = 4
            ax.add_patch(patches.Rectangle((j + 0.5, i + 0.1), 0.9, 0.9, facecolor=facecolor, edgecolor=edgecolor))
            ax.text(j + 0.6, i + 0.4, label, fontsize=fontsize, color="white")

    ax.set_xlim(0, len(seq._seq) + 1)
    ax.set_ylim(0, len(msa) + 0.1)

    custom_xticks = [i + 1 for i in range(len(seq._seq))]
    ax.set_xticks(custom_xticks)
    ax.set_xticklabels(custom_xticks)

    custom_yticks = [i + 0.6 for i in range(len(msa))]
    ax.set_yticks(custom_yticks)
    ax.set_yticklabels(seq_labels)
    ax.tick_params(axis=u'both', which=u'both', length=0)

    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    ax.set_aspect("equal")
    plt.subplots_adjust(left=0.3)
    plt.savefig(path, dpi=300)
    plt.clf()
    plt.close()

def main() -> None:
    """
    Main function.
    """
    args = cli()

    # Parse sequences.
    names, primary_sequences, binned_sequences = [], [], []
    with open(args.input, "r") as file_open:
        for i, line in enumerate(file_open):
            # Every first line with '>' is header, second line is primary sequence, and third line is binned sequence.
            if line.startswith(">"):
                names.append(line.strip()[1:].split("|")[0])
                primary_sequences.append(next(file_open).strip().split("|"))
                binned_sequences.append(next(file_open).strip().split("|"))
            print(f"{i}".zfill(6), end="\r")

    # Create fingerprints.
    X = []
    for i, seq in enumerate(binned_sequences):
        X.append(get_biosynthetic_fingerprint(seq))
        print(f"{i}".zfill(6), end="\r")
    X = np.array(X)
    print("Shape of fingerprints:", X.shape)

    # Perform UMAP on fingerprints.
    reducer = umap.UMAP(n_neighbors=10, min_dist=0.1, metric="cosine", random_state=42)
    embedding = reducer.fit_transform(X)
    print("Shape of embedding:", embedding.shape)

    # Cluster original data with HDBSCAN.
    model = HDBSCAN(min_cluster_size=10, min_samples=10)
    model.fit(embedding)
    num_clusters = len(set(model.labels_))

    # Visualie umap in 2D.
    cmap = plt.cm.get_cmap("tab20")    
    plt.scatter(embedding[:, 0], embedding[:, 1], s=0.1, c=model.labels_, cmap=cmap)
    plt.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
    plt.title(f"{X.shape[0]} sequences / {num_clusters} clusters")
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.gca().set_facecolor("black")
    plt.savefig(f"{args.output}/umap_biosynthetic_space.png", dpi=300)
    plt.clf()

    # Save embedding.
    np.save(f"{args.output}/embedding.npy", embedding)

    # Save names.
    with open(f"{args.output}/embedding_identifiers.txt", "w") as file_open:
        for name in names:
            file_open.write(f"{name}\n")

    # Save embedder.
    joblib.dump(reducer, f"{args.output}/embedder.joblib")

    # Per cluster, get names, primary sequences, and binned sequences.
    for i in sorted(list(set(model.labels_))):
        if i > 0:
            # Get indices of cluster.
            indices = np.where(model.labels_ == i)[0]
            print(f"Cluster {i} has {len(indices)} sequences.")

            # Get names, primary sequences, and binned sequences.
            names_cluster = [names[index] for index in indices]
            primary_sequences_cluster = [primary_sequences[index] for index in indices]
            binned_sequences_cluster = [binned_sequences[index] for index in indices]

            # Align and save alignment figure.
            seqs = [ModuleSequence(name, binned_seq) for name, binned_seq in zip(names_cluster, binned_sequences_cluster)]
            msa = MultipleSequenceAlignment(seqs, gap_cost=5, gap_end_cost=2)

            path = f"{args.output}/cluster_{i}"
            plot_msa(path, msa.msa[list(msa.msa.keys())[0]])

        print(f"{i}".zfill(6), end="\r")

    exit(0)

if __name__ == "__main__":
    main()