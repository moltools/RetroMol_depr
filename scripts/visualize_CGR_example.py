#!/usr/bin/env python3
"""
Description:    This script visualizes the CGR space.
Usage:          python3 visualize_CGR_example.py -o <path/to/output/dir>
Dependencies:   imageio (pip install imageio, see for more info: https://stackoverflow.com/questions/753190/programmatically-generate-video-or-animated-gif-in-python)
"""
import argparse
import os

import imageio
import matplotlib.pyplot as plt

from retromol_sequencing.fingerprint import CompoundClassMapping

def cli() -> argparse.Namespace:
    """
    Parse command line arguments.

    :return: Arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(description="Visualize CGR example.")
    parser.add_argument("-o", "--output", type=str, help="Output directory.", required=True)
    return parser.parse_args()

def get_CGR_space(show_labels: bool) -> None:
    """
    Get CGR space.

    :param bool show_labels: Show labels.
    :return: Figure.
    :rtype: plt.Figure
    """
    name_to_label = {
        "A": "A",
        "B": "B",
        "C": "C",
        "D": "D",
        "PolarAndCharged": "Polar and charged",
        "SmallHydrophobic": "Small hydrophobic",
        "SmallNonHydrophobic": "Small non-hydrophobic",
        "Tiny": "Tiny",
        "Bulky": "Bulky (mainly phenyl derivatives)",
        "CyclicAliphatic": "Cyclic aliphatic",  
        "SulfurContaining": "Cysteine-like",
        "Undefined": "Undefined",
    }

    circle = plt.Circle((0, 0), 1, fill=False, linestyle="--", alpha=0.25)    
    plt.gca().add_patch(circle)

    for motif, coords in CompoundClassMapping.__members__.items():
        plt.scatter(*coords.value, c="k", marker="o")

        if show_labels:
            x, y = coords.value
            plt.text(x + 0.025, y + 0.025, name_to_label[motif], fontsize=8)
            

    plt.gca().set_aspect("equal", adjustable="box")
    plt.xticks([])
    plt.yticks([])
    plt.axis("off")    

def main() -> None:
    """
    Main function.
    """
    args = cli()

    get_CGR_space(show_labels=True)
    plt.savefig(os.path.join(args.output, "step_000.png"), dpi=300)
    plt.clf()

    # Define example sequence.
    seq = ["B", "B", "A", "D", "B", "D", "C", "C", "Bulky"]

    locs = [(0, 0)]

    prev_loc = (0, 0)
    for i, motif in enumerate(seq):
        next_loc = CompoundClassMapping.get_mapping(motif)        
        middle_loc = ((prev_loc[0] + next_loc[0]) / 2, (prev_loc[1] + next_loc[1]) / 2)
        locs.append(middle_loc)
        prev_loc = middle_loc

        get_CGR_space(show_labels=False)
        plt.plot([x[0] for x in locs], [x[1] for x in locs], c="r", marker="o", markersize=2)
        plt.text(-1.5, -1.0, f"{motif}", bbox=dict(facecolor="yellow", alpha=0.5, edgecolor="black"))
        plt.savefig(os.path.join(args.output, "step_" + f"{i + 1}".zfill(3) + ".png"), dpi=300)
        plt.clf()

    # Read and sort all image png paths from out dir. 
    paths = sorted([os.path.join(args.output, x) for x in os.listdir(args.output) if x.endswith(".png")])
    images = []
    for filename in paths:
        images.append(imageio.imread(filename))
    imageio.mimsave(os.path.join(args.output, "CGR_example.gif"), images, fps=1)
        
    exit(0)

if __name__ == "__main__":
    main()