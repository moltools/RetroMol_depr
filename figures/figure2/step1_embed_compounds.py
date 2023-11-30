#!/usr/bin/env python3
"""
Embed compounds in 2D space using UMAP.

Dependencies: 
    - matplotlib
    - numpy
    - rdkit
    - tqdm
    - umap-learn
"""
import argparse 
import matplotlib.pyplot as plt
import numpy as np
import os
import umap 
from tqdm import tqdm 
from rdkit import Chem
from rdkit import DataStructs

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, required=True, help="Path to compounds TSV file.")
parser.add_argument("--header", action="store_true", help="Input file contains a header line.")
args = parser.parse_args()

# Savefile paths.
lbl_path = "figures/figure2/labels.npy"
fps_path = "figures/figure2/fps.npy"
emb_path = "figures/figure2/embedding.npy"

if not os.path.exists(lbl_path) and not os.path.exists(fps_path):
    # Parse compounds from input file.
    labels, compounds = [], []
    with open(args.input, "r") as fo:
        if args.header: next(fo)
        for line in tqdm(fo, desc="Parsing records from input file"):
            line = line.strip()
            npaid, smiles = line.split("\t")
            labels.append(npaid)
            compounds.append(smiles)

    # Convert compounds to molecular fingerprints.
    fps = []
    for compound in tqdm(compounds, desc="Converting compounds to molecular fingerprints"):
        mol = Chem.MolFromSmiles(compound)
        fp = Chem.RDKFingerprint(mol, maxPath=5, fpSize=2048)
        fps.append(np.array(fp))
    X = np.array(fps)
    print("Shape of fingerprint data:", X.shape)

    # Save data.
    print("Saving data...")
    np.save(lbl_path, labels)
    np.save(fps_path, X)

else:
    # Load saved data.
    print("Loading saved data...")
    labels = np.load(lbl_path)
    X = np.load(fps_path)
    print("Shape of fingerprint data:", X.shape)

# Create UMAP embedding.
print("Creating UMAP embedding...")
embedder = umap.UMAP(n_neighbors=50, metric="cosine")
embedding = embedder.fit_transform(X)
print("Embedding shape:", embedding.shape)

# View embedding.
plt.scatter(embedding[:, 0], embedding[:, 1], s=1)
plt.show()

# Save embedding.
print("Saving embedding...")
np.save(emb_path, embedding)

exit("Done.")