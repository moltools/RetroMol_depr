#!/usr/bin/env python3
"""
Description:    This script trains a classifier for the class predictions of motifs for sequence alignment.
Usage:          python3 train_classifier.py -o <path/to/output/dir>
"""
import argparse
import joblib
import os 

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from rdkit import Chem
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import StandardScaler

from retromol_sequencing.fingerprint import get_amino_acid_fingerprint

classes = [
    "polar and charged",
    "small hydrophobic",
    "bulky (mainly phenyl derivatives)",
    "small non-hydrophobic",
    "cyclic aliphatic",
    "tiny",
    "cysteine-like",
]

# Source: https://srobinson.shinyapps.io/AdenylPred/#section-substrate-groups
compounds_with_classes = [
    (r"C(CC(C(=O)O)N)CN=C(N)N", "polar and charged"),
    (r"C(CC(=O)N)C(C(=O)O)N", "polar and charged"),
    (r"C(CC(=O)[O])C(C(=O)[O])[NH2]", "polar and charged"),
    (r"C(C(C(=O)[O])N)C(=O)[O]", "polar and charged"),
    (r"C(C(C(=O)O)N)O", "polar and charged"),
    (r"CC(C)C(C(=O)O)N", "small hydrophobic"),
    (r"CCC(C)C(C(=O)O)N", "small hydrophobic"),
    (r"CC(C)C(C(=O)O)O", "small hydrophobic"),
    (r"CCC(C(=O)O)N", "small hydrophobic"),
    (r"CC(C)C(CO)N", "small hydrophobic"),
    (r"CC1CC1(C(=O)O)N", "small hydrophobic"),
    (r"CC(C)CC(=O)C(=O)[O]", "small hydrophobic"),
    (r"CC(C)C(=O)C(=O)[O]", "small hydrophobic"),
    (r"C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N", "bulky (mainly phenyl derivatives)"),
    (r"C1=CC=C(C=C1)CC(C(=O)O)N", "bulky (mainly phenyl derivatives)"),
    (r"C1=CC(=CC=C1CC(C(=O)O)N)O", "bulky (mainly phenyl derivatives)"),
    (r"C(CCN)CC(C(=O)O)N", "bulky (mainly phenyl derivatives)"),
    (r"C1=C(NC=N1)CC(C(=O)O)N", "bulky (mainly phenyl derivatives)"),
    (r"C1=CC=C(C(=C1)C(=O)CC(C(=O)O)N)N", "bulky (mainly phenyl derivatives)"),
    (r"CC(C(=O)O)N", "small non-hydrophobic"),
    (r"CC(C(C(=O)O)N)O", "small non-hydrophobic"),
    (r"C1CC(NC1)C(=O)O", "cyclic aliphatic"),
    (r"C1CC[NH]C(C1)C(=O)[O]", "cyclic aliphatic"),
    (r"CC1CNC(C1)C(=O)O", "cyclic aliphatic"),
    (r"CC=CC1CC(NC1)C(=O)O", "cyclic aliphatic"),
    (r"C1C(CNC1C(=O)O)O", "cyclic aliphatic"),
    (r"CCCC1CC(NC1)C(=O)O", "cyclic aliphatic"),
    (r"CC1(CCNC1C(=O)O)O", "cyclic aliphatic"),
    (r"C(C(=O)O)N", "tiny"),
    (r"C(C(=O)[O])O", "tiny"),
    (r"C(C(C(=O)O)N)S", "cysteine-like"),
]

def cli() -> argparse.Namespace:
    """
    Command line interface.

    :returns: Command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", type=str, required=True, help="Output directory.")
    return parser.parse_args()

def main() -> None:
    """
    Main function.
    """
    args = cli()

    # Prepare training data.
    X, y = [], []
    for smiles, assigned_class in compounds_with_classes:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        X.append(get_amino_acid_fingerprint(mol))
        y.append(classes.index(assigned_class))
    X, y = np.array(X), np.array(y)
    print(X.shape, y.shape)

    # Visualize training data.
    num_components = 2
    pca = PCA(n_components=num_components)
    scaler = StandardScaler()
    X_norm = scaler.fit_transform(X)
    pcs = pca.fit_transform(X_norm)

    plt.figure(figsize=(8, 4))
    colors = ["b", "g", "r", "c", "m", "y", "k"]
    legend_labels = [classification.capitalize() for classification in classes]
    for i in range(len(classes)):
        plt.scatter(pcs[y == i][:, 0], pcs[y == i][:, 1], c=colors[i], marker="o", label=legend_labels[i])
    legend_patches = [mpatches.Patch(color=colors[i], label=legend_labels[i]) for i in range(len(classes))]
    plt.legend(handles=legend_patches, bbox_to_anchor=(1.02, 1.0))
    ev = pca.explained_variance_ratio_
    plt.xlabel(f"PC 1 ({round(ev[0] * 100, 2)}%)")
    plt.ylabel(f"PC 2 ({round(ev[1] * 100, 2)}%)")
    plt.tight_layout()
    out_path = os.path.join(args.output, "training_data.png")
    plt.savefig(out_path, dpi=300)
    plt.clf()
    
    # Train classifier.
    model = RandomForestClassifier(n_estimators=1000, random_state=42)
    model.fit(X, y)
    y_pred = model.predict(X)
    accuracy = accuracy_score(y, y_pred)
    print(f"Accuracy: {accuracy}")
    out_path = os.path.join(args.output, "model.joblib")
    joblib.dump(model, out_path)

    # Visualize feature importance in bar plot.
    feat_labels = [
        "Number of heavy atoms",
        "Number of sulfur atoms",
        "Number of 3-atom-sized rings",
        "Number of 5-atom-sized rings",
        "Number of 6-atom-sized rings",
        "Number of hydrogen bond acceptors",
        "Number of hydrogen bond donors",
        "Bertz graph complexity index",
        "Number of N-H and O-H bonds",
        "Molecular weight",
        "LogP",
    ]
    importances = model.feature_importances_
    indices = np.argsort(importances)[::-1]
    plt.figure(figsize=(4, 8))
    plt.bar(range(X.shape[1]), importances[indices], color="k", align="center")
    plt.xticks(range(X.shape[1]), [feat_labels[i] for i in indices], rotation=90)
    plt.xlim([-1, X.shape[1]])
    plt.tight_layout()
    out_path = os.path.join(args.output, "feature_importances.png")
    plt.savefig(out_path, dpi=300)
    plt.clf()

    exit(0)

if __name__ == "__main__":
    main()