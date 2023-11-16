#!/usr/bin/env python
import argparse
import joblib
import os

import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem, RDLogger
from sklearn.manifold import TSNE 
from sklearn.decomposition import PCA

from retromol.chem import mol_to_fingerprint

def cli() -> argparse.Namespace:
    """
    Command-line interface.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="Path to dataset file with compounds.")
    parser.add_argument("-m", "--model", type=str, required=True, help="Path to model save file.")
    parser.add_argument("-o", "--output-dir", type=str, required=True, help="Path to output dir.")
    return parser.parse_args()

def main() -> None:
    """
    Driver code.
    """
    RDLogger.DisableLog("rdApp.*")
    args = cli()

    seqs = []
    smis = []
    with open(args.input, "r") as f:
        for line in f:
            smi, seq = line.strip().split(",")
            smis.append(smi)
            seqs.append(seq)

        # # select first 1000
        # smis = smis[:1000]
    
    with open(args.model, "rb") as f:
        model = joblib.load(f)

    X = []
    has_cycle = []
    for i, smi in enumerate(smis):
        mol = Chem.MolFromSmiles(smi)
        if mol.GetRingInfo().NumRings() > 0:
            has_cycle.append(1)
        else:
            has_cycle.append(0)
        fp = mol_to_fingerprint(mol, 4, 2048)
        X.append(fp)
        print(f"{i}".zfill(5), end="\r")
    X = np.array(X)
    has_cycle = np.array(has_cycle)

    print(X.shape)

    preds = model.predict_proba(X)[:, 1]

    # Get top 10 predictions.
    print("max pred:", np.max(preds))
    idx = np.argsort(preds)[::-1][:10]
    for i in idx:
        print(smis[i], preds[i])

    for i in idx:
        print(seqs[i])

    # Get bottom 10 predictions.
    print("min pred:", np.min(preds))
    idx = np.argsort(preds)[:10]
    for i in idx:
        print(smis[i], preds[i])

    # # PCA embedding.
    # pca = PCA(n_components=2, random_state=42)
    # X_pca = pca.fit_transform(X)
    # imps = pca.explained_variance_ratio_

    # sc = plt.scatter(X_pca[:, 0], X_pca[:, 1], c=preds, cmap="coolwarm", alpha=0.5)
    # cbar = plt.colorbar(sc)
    # cbar.set_label("Predicted to be antibacterial [0, 1]")
    # plt.xlabel("PCA dimension 1 (EV={:.2f})".format(imps[0]))
    # plt.ylabel("PCA dimension 2 (EV={:.2f})".format(imps[1]))
    # path = os.path.join(args.output_dir, "pca_1.png")
    # plt.savefig(path, dpi=300)
    # plt.clf()

    # # Define two colors and plot if the molecule has a cycle or not.
    # X_has_cycle = X_pca[has_cycle == 1]
    # X_no_cycle = X_pca[has_cycle == 0]

    # sc = plt.scatter(X_has_cycle[:, 0], X_has_cycle[:, 1], c="green", alpha=0.5, label="Has cycle")
    # sc = plt.scatter(X_no_cycle[:, 0], X_no_cycle[:, 1], c="red", alpha=0.5, label="No cycle")

    # plt.xlabel("PCA dimension 1 (EV={:.2f})".format(imps[0]))
    # plt.ylabel("PCA dimension 2 (EV={:.2f})".format(imps[1]))
    # path = os.path.join(args.output_dir, "pca_2.png")
    # plt.savefig(path, dpi=300)

    # # TSNE embedding.   
    # tsne = TSNE(n_components=2, random_state=42)
    # X_tsne = tsne.fit_transform(X)

    # sc = plt.scatter(X_tsne[:, 0], X_tsne[:, 1], c=preds, cmap="coolwarm", alpha=0.5)
    # cbar = plt.colorbar(sc)
    # cbar.set_label("Predicted to be antibacterial [0, 1]")
    # plt.xlabel("t-SNE dimension 1")
    # plt.ylabel("t-SNE dimension 2")
    # path = os.path.join(args.output_dir, "tsne_1.png")
    # plt.savefig(path, dpi=300)
    # plt.clf()

    # # Define two colors and plot if the molecule has a cycle or not.
    # X_has_cycle = X_tsne[has_cycle == 1]
    # X_no_cycle = X_tsne[has_cycle == 0]

    # sc = plt.scatter(X_has_cycle[:, 0], X_has_cycle[:, 1], c="green", alpha=0.5, label="Has cycle")
    # sc = plt.scatter(X_no_cycle[:, 0], X_no_cycle[:, 1], c="red", alpha=0.5, label="No cycle")

    # plt.xlabel("t-SNE dimension 1")
    # plt.ylabel("t-SNE dimension 2")
    # path = os.path.join(args.output_dir, "tsne_2.png")
    # plt.savefig(path, dpi=300)

    exit(0)

if __name__ == "__main__":
    main()
