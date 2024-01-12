#!/usr/bin/env python3
import argparse
import math
import os
import typing as ty
from dataclasses import dataclass

import numpy as np
import umap
import plotly.express as px 

def cli() -> argparse.Namespace:
    """
    Command line interface.
    
    :return: argparse.Namespace
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--fasta", type=str, required=True, help="Path to FASTA file containing primary sequences.")
    parser.add_argument("-o", "--out", type=str, required=True, help="Path to output dir.")
    return parser.parse_args()

@dataclass
class Record:
    """
    Dataclass for FASTA record.

    :param name: ID of record.
    :type name: str
    :param header: Header of record.
    :type header: str
    """
    name: str
    seq: ty.List[str]

def read_fasta(fasta: str) -> ty.List[Record]:
    """
    Read FASTA file.
    
    :param fasta: Path to FASTA file.
    :type fasta: str
    :return: List of FASTA records.
    :rtype: ty.List[Record]
    """
    records = []

    with open(fasta, "r") as file_open:
        for line in file_open:
            line = line.strip()

            # Parse header.
            if line.startswith(">"):
                name = line[1:].split("|")[0]

            # Parse sequence.
            else:
                seq = line.split("|")
                records.append(Record(name, seq))

    return records

def sequence_to_class(seq: ty.List[str]) -> str:
    assigned = []
    for x in seq:
        if x[0].isalpha() and x[1:].isdigit() and x[0] == "A": assigned.append("PK")
        if x[0].isalpha() and x[1:].isdigit() and x[0] == "B": assigned.append("PK")
        if x[0].isalpha() and x[1:].isdigit() and x[0] == "C": assigned.append("PK")
        if x[0].isalpha() and x[1:].isdigit() and x[0] == "D": assigned.append("PK")
        else: assigned.append("NRP")
    if all([x == "PK" for x in assigned]): return "PK"
    if all([x == "NRP" for x in assigned]): return "NRP"
    else: return "Hybrid"

def compound_to_position(compound: str) -> ty.Tuple[float, float]:
    """
    Convert compound to position.
    
    :param compound: Compound.
    :type compound: str
    :return: Position.
    :rtype: ty.Tuple[float, float]
    """
    mapping = {
        "A": (1.000, 0.000), 
        "B": (0.841, 0.541),
        "C": (0.415, 0.910), 
        "D": (-0.142, 0.990),
        "polar_and_charged": (-0.655, 0.756),
        "small_hydrophobic": (-0.959, 0.282),
        "bulky_mainly_phenyl_derivatives": (-0.959, -0.282),
        "small_non_hydrophobic": (-0.655, -0.756),
        "cyclic_aliphatic": (-0.142, -0.990),
        "cysteine": (0.415, -0.910),
        "tiny": (0.841, -0.541),
    }
    x = compound

    if x[0].isalpha() and x[1:].isdigit() and x[0] == "A": return mapping["A"]
    if x[0].isalpha() and x[1:].isdigit() and x[0] == "B": return mapping["B"]
    if x[0].isalpha() and x[1:].isdigit() and x[0] == "C": return mapping["C"]
    if x[0].isalpha() and x[1:].isdigit() and x[0] == "D": return mapping["D"]
    if x in ["Arg", "Gln", "Aad", "Asp", "aAsp/Asp", "Glu", "aGlu/Glu", "Asn", "OH-Asn", "MeO-Glu", "bOMe-Asp", "Dab", "Ser", "Ph-Ser", "Ac-Ser"] or "Orn" in x: return mapping["polar_and_charged"]
    if x in ["Val", "bOH-Val", "Aib", "Leu", "3OH-Leu", "Ile/aIle", "NFo-Leu", "NAc-Leu"]: return mapping["small_hydrophobic"] 
    if x in ["Tyr", "bOH-Tyr", "Phe", "Trp", "OH-Trp", "dh-Trp", "Hpg"]: return mapping["bulky_mainly_phenyl_derivatives"]
    if x in ["aThr/Thr", "Ala", "Abu", "dhAbu", "bAla", "Dpr"]: return mapping["small_non_hydrophobic"]
    if x in ["Pro", "NFo-Pro", "aPro/Pro", "NMe-Hpr", "Hpr"]: return mapping["cyclic_aliphatic"]
    if x in ["Cys", "NAc-Cys", "Met"]: return mapping["cysteine"]
    if x in ["Gly"]: return mapping["tiny"]
    else: 
        raise ValueError(f"Invalid compound {compound}")

def create_fingerprint(seq: ty.List[str]) -> np.ndarray:
    """
    Create fingerprint for sequence.
    
    :param seq: Sequence.
    :type seq: ty.List[str]
    :return: Fingerprint.
    :rtype: np.ndarray

    See also: https://srobinson.shinyapps.io/AdenylPred/#section-substrate-groups
    """
    fp = [0 for _ in range(3600)]

    def find_bin_index(x, y):
        # Get clockwise angle from x-axis.
        angle = math.degrees(math.atan2(y, x))
        if angle < 0: angle += 360
        angle = math.ceil(angle)

        # Get distance to origin.
        dist = math.ceil(math.sqrt(x**2 + y**2) * 10)

        # Get bin index.
        bin_index = angle * dist 

        return bin_index
        
    prev_loc = (0, 0)
    for motif in seq:
        next_loc = compound_to_position(motif)
        middle_loc = ((prev_loc[0] + next_loc[0]) / 2, (prev_loc[1] + next_loc[1]) / 2)
        bin_index = find_bin_index(middle_loc[0], middle_loc[1])
        fp[bin_index] += 1
        prev_loc = middle_loc

    fp = np.array(fp)
    fp = fp / np.linalg.norm(fp)

    return fp

def main() -> None:
    """
    Main function.
    """
    # Parse command line arguments.
    args = cli()
    records = read_fasta(args.fasta)
    print(f"Number of records: {len(records)}")

    # Create fingerprints.
    parsed_records = []
    labels, fps, classes = [], [], []
    for i, record in enumerate(records):
        try:
            fp = create_fingerprint(record.seq)
            fps.append(fp)
            labels.append(record.name)
            classes.append(sequence_to_class(record.seq))
            parsed_records.append(record)

        except Exception:
            pass

        padding = len(str(len(records)))
        print(f"{i}".zfill(padding), end="\r")

    # To numpy arrays.
    fps = np.array(fps)

    # Print dimensionality.
    print(f"Dimensionality: {fps.shape[0]} by {fps.shape[1]} (num_labels={len(labels)})")

    # UMAP.
    reducer = umap.UMAP(n_neighbors=15, min_dist=1.0, metric="euclidean")
    embedding = reducer.fit_transform(fps)
    print(f"Embedding shape: {embedding.shape}")

    # Give colors to selected labels.
    colors = []
    for label in labels:
        if not label.startswith("NPA"):
            colors.append("Selected")
        else:
            colors.append("Other")

    # Plot with Plotly.
    fig = px.scatter(
        embedding, 
        x=embedding[:, 0], 
        y=embedding[:, 1], 
        hover_name=labels,
        # color=colors,
        # color=classes
    )
    fig.write_html(os.path.join(args.out, "embedding.html"))

    # Save fingerprints.
    np.save(os.path.join(args.out, "fingerprints.npy"), fps)
    # Save labels.
    with open(os.path.join(args.out, "labels.txt"), "w") as file_open:
        for label in labels:
            file_open.write(f"{label}\n")

    # Save embedding.
    np.save(os.path.join(args.out, "embedding.npy"), embedding)

    # Write out new fasta with only selected sequences in same order as fingerprints.
    with open(os.path.join(args.out, "selected.fasta"), "w") as file_open:
        for i, record in enumerate(parsed_records):
            file_open.write(f">{record.name}\n")
            file_open.write(f"{'|'.join(record.seq)}\n")

    exit(0)

if __name__ == "__main__":
    main()