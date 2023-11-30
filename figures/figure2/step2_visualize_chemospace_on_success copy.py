"""
Visualize the chemospace of parsed compounds.
"""
import argparse
import glob 
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm 

from retromol.parsing import Result 

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, required=True, help="Path to dir with RetroMol JSON results.")
args = parser.parse_args()

# File paths.
emb_path = "figures/figure2/embedding.npy"
lbl_path = "figures/figure2/labels.npy"
npatlas_path = "npatlas/NPAtlas_download_20231129.json"

# Load data.
X = np.load(emb_path)
print("Shape of embedding data:", X.shape)
labels_names = np.load(lbl_path)
print("Shape of labels data:", labels_names.shape)

# Read RetroMol results JSONs.
filepaths = glob.glob(args.input + "/*.json")

# Parse results.
# Get empty array of success labels of length of labels_names.
labels_success = np.zeros(len(labels_names))

for filepath in tqdm(filepaths, desc="Reading RetroMol results"):
    result = Result.from_json(filepath)
    if result.name not in labels_names: continue
    else: 
        # Get index of result in labels_names.
        index = np.where(labels_names == result.name)[0][0]
        if result.success:
            total_atoms = result.substrate.GetNumAtoms()
            score = result.score
            score_as_percentage = round(score / total_atoms * 100)
            if score == 0:
                labels_success[index] = 1
            elif score_as_percentage == 100:
                labels_success[index] = 2
            else:
                labels_success[index] = 3
        else:
            labels_success[index] = 0
labels_success = np.array(labels_success)
print("Shape of labels_success data:", labels_success.shape)

# Plot embedding with success labels.
plt.scatter(X[labels_success == 2, 0], X[labels_success == 2, 1], c="orange", alpha=0.75, s=5, label="Fully unidentified")
plt.scatter(X[labels_success == 3, 0], X[labels_success == 3, 1], c="blue", alpha=0.75, s=5, label="Partially unidentified")
plt.scatter(X[labels_success == 0, 0], X[labels_success == 0, 1], c="red", alpha=0.75, s=5, label="Failed")
plt.scatter(X[labels_success == 1, 0], X[labels_success == 1, 1], c="green", alpha=0.75, s=5, label="Fully identified")
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
plt.savefig("figures/figure2/chemospace_success.png", dpi=300, bbox_inches="tight")