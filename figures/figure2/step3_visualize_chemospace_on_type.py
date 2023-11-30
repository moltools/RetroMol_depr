"""
Visualize the chemospace of parsed compounds.
"""
import glob 
import json
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm 

# File paths.
emb_path = "figures/figure2/embedding.npy"
lbl_path = "figures/figure2/labels.npy"
npatlas_path = "npatlas/NPAtlas_download_20231129.json"

# Load data.
X = np.load(emb_path)
print("Shape of embedding data:", X.shape)
labels_names = np.load(lbl_path)
print("Shape of labels data:", labels_names.shape)

# Parse results.
# Get empty array of success labels of length of labels_names.
labels_success = np.zeros(len(labels_names))

data = json.load(open(npatlas_path))
print("Number of compounds in dataset:", len(data)) # Expected: 33372 compounds
for record in tqdm(data, desc="Parsing records from NPAtlas"):
    try:
        npaid = record["npaid"]
        smiles = record["smiles"]
        pathway_classification = record["npclassifier"]["pathway_results"]

        if npaid in labels_names:
            if "Polyketides" in pathway_classification and "Amino Acids and Peptides" in pathway_classification:
                labels_success[np.where(labels_names == npaid)[0]] = 0
            elif "Polyketides" in pathway_classification:
                labels_success[np.where(labels_names == npaid)[0]] = 1
            elif "Amino acids and Peptides" in pathway_classification:
                labels_success[np.where(labels_names == npaid)[0]] = 2
            else:
                raise ValueError("Compound not in allowed pathway classifications.")

    except Exception:
        continue
print("Labels success shape:", labels_success.shape)

# Plot embedding with success labels.
plt.scatter(X[labels_success == 1, 0], X[labels_success == 1, 1], c="green", alpha=0.75, s=5, linewidths=0, label="Amini Acids and Peptides")
plt.scatter(X[labels_success == 2, 0], X[labels_success == 2, 1], c="blue", alpha=0.75, s=5, linewidths=0, label="Both")
plt.scatter(X[labels_success == 0, 0], X[labels_success == 0, 1], c="red", alpha=0.75, s=5, linewidths=0, label="Polyketides")
plt.title(f"UMAP of compound fingerprints (n={len(labels_success)})")
plt.xlabel("UMAP component 1")
plt.ylabel("UMAP component 2")
plt.xticks([])
plt.yticks([])
plt.gca().axes.xaxis.set_ticklabels([])
plt.gca().axes.yaxis.set_ticklabels([])
plt.gca().axes.xaxis.set_ticks([])
plt.gca().axes.yaxis.set_ticks([])
plt.legend() 
plt.savefig("figures/figure2/chemospace_type.png", dpi=300, bbox_inches="tight")