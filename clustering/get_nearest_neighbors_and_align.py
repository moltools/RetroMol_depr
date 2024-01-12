#!/usr/bin/env python3
import argparse 
import typing as ty
import re
import os

import numpy as np
from sklearn.neighbors import KDTree

import matplotlib.pyplot as plt
import matplotlib.patches as patches

from cluster_sequences import create_fingerprint
from aligner.parser import Record 
from aligner.aligner import multiple_sequence_alignment
from aligner.aligner import PKModule

def cli() -> argparse.Namespace:
    """
    Command line interface.
    
    :return: argparse.Namespace
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--fingerprints", type=str, required=True, help="Path to fingerprints file.")
    parser.add_argument("-l", "--labels", type=str, required=True, help="Path to labels file.")
    parser.add_argument("-s", "--sequences", type=str, required=True, help="Corresponding sequences in fasta format.")
    parser.add_argument("-o", "--out", type=str, required=True, help="Path to output dir.")
    return parser.parse_args()

def get_nearest_neighbors(fingerprints: np.ndarray, query: np.ndarray, k: int) -> ty.Tuple[np.ndarray, np.ndarray]: 
    """
    Use KD tree to get nearest neighbors.

    :param fingerprints: Fingerprints.
    :type fingerprints: np.ndarray
    :param query: Query.
    :type query: np.ndarray
    :param k: Number of nearest neighbors.
    :type k: int
    :return: Nearest neighbors with distances.
    :rtype: ty.Tuple[np.ndarray, np.ndarray]
    """
    tree = KDTree(fingerprints, leaf_size=2)
    dist, ind = tree.query(query, k=k)
    return dist[0], ind[0]

def main() -> None:
    """
    Main function.
    """
    args = cli()

    # Read in fingerprints. Its a npy file.
    fingerprints = np.load(args.fingerprints)
    print(fingerprints.shape)

    # Read labels. Label per line.
    labels = []
    with open(args.labels, "r") as file_open:
        for line in file_open:
            line = line.strip()
            labels.append(line)
    labels = np.array(labels)
    print(labels.shape)

    # Define query.
    # query = ["C1", "C1", "B2", "B1", "D2", "D2", "B2", "C1", "B1", "B2", "C1", "C1"]
    query = ["D2", "B1", "A2", "C2", "B5", "A2", "D2", "C1", "C1", "C2", "B1", "A1", "D2", "B11"]
    fp = np.array([create_fingerprint(query)])

    # Get nearest neighbors.
    dists, nbs = get_nearest_neighbors(fingerprints, fp, 10)
    
    for i, nb in enumerate(nbs):
        dist = dists[i]
        label = labels[nb]
        print(f">> ({i}) {label} ({round(dist, 3)})")

    # Get corresponding sequences from nearest neighbors.
    records = []
    with open(args.sequences, "r") as file_open:
        for line in file_open:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split("|")[0]
            else:
                seq = line.split("|")
                records.append(Record(name, seq))
    print(len(records))

    nbs_seqs = []
    for i, nb in enumerate(nbs):
        nbs_seqs.append(records[nb])

    # Add query.
    nbs_seqs.append(Record("Query", query))
    dists = np.append(dists, 0)

    # Now we have nbs, dists, and seqs.
    print(len(dists), len(nbs), len(nbs_seqs))

    lookup = {record.name: dists[i] for i, record in enumerate(nbs_seqs)}

    # Align sequences.
    msa = multiple_sequence_alignment(nbs_seqs, gap_cost=5, gap_end_cost=1)

    # msa.display()

    # aligned_seqs = msa.msa[list(msa.msa.keys())[0]]
    # aligned_names = [aligned_seq.name for aligned_seq in aligned_seqs]
    # aligned_original_seqs = [(name, lookup[name][0], lookup[name][1], aligned_seqs[i]) for i, name in enumerate(aligned_names)]

    # TODO: put query in list for alignment

    labels = []
    dist_labels = []
    aligned_seqs = []
    for seq in msa.msa[list(msa.msa.keys())[0]]:
        label = f"{seq.name}"
        labels.append(label)
        dist_label = f"{round(lookup[seq.name], 3)}"
        dist_labels.append(dist_label)
        new_seq = []
        for motif in seq._seq:
            motif = motif[0]
            if motif == PKModule.GAP: new_seq.append("GAP")
            else: new_seq.append(motif.original_name)
        aligned_seqs.append(new_seq)
        print(label, dist_label, new_seq)


    # scale figure size with number and length of sequences
    figsize_x = len(aligned_seqs[0])
    figsize_y = len(aligned_seqs)

    fig, ax = plt.subplots(figsize=(figsize_x, figsize_y))

    for i, seq in enumerate(aligned_seqs):
        for j, item in enumerate(seq):
            if item == "GAP":
                facecolor, edgecolor = "white", "white"
                label = ""
            elif re.match(r"A\d", item) or (re.match(r"A", item) and len(item) == 1):
                facecolor, edgecolor = "red", "black"
                label = item
            elif re.match(r"B\d", item) or re.match(r"B", item):
                facecolor, edgecolor = "blue", "black"
                label = item
            elif re.match(r"C\d", item) or re.match(r"C", item):
                facecolor, edgecolor = "green", "black"
                label = item
            elif re.match(r"D\d", item) or re.match(r"D", item):
                facecolor, edgecolor = "orange", "black"
                label = item
            else:
                label = item
                facecolor, edgecolor = "grey", "black"
                label = ""
            ax.add_patch(patches.Rectangle((j + 0.5, i + 0.1), 0.9, 0.9, facecolor=facecolor, edgecolor=edgecolor))
            ax.text(j + 0.6, i + 0.4, label, fontsize=12, color="white")

    ax.set_xlim(0, len(seq) + 1)
    ax.set_ylim(0, len(aligned_seqs) + 0.1)
    # ax.set_xlabel("Position")
    # ax.set_ylabel("Sequence")
    # ax.set_title("MSA of community")

    custom_xticks = [i + 1 for i in range(len(seq))]
    ax.set_xticks(custom_xticks)
    ax.set_xticklabels(custom_xticks)

    custom_yticks = [i + 0.6 for i in range(len(aligned_seqs))]
    ax.set_yticks(custom_yticks)
    ax.set_yticklabels(labels)

    # Also put labels on right side of plot y axis
    # ax2 = ax.twinx()
    # ax2.set_yticks([i for i in range(len(aligned_seqs))])
    # ax2.set_yticklabels(dist_labels)
    # ax2.tick_params(axis=u'both', which=u'both', length=0)

    ax.tick_params(axis=u'both', which=u'both', length=0)

    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    ax.set_aspect("equal")

    plt.subplots_adjust(left=0.3)
    
    # plt.show()
    plt.savefig(os.path.join(args.out, "alignment.png"), dpi=300)

    exit(0) 

    exit(0)

if __name__ == "__main__":
    main()