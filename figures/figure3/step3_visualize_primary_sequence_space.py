"""
Visualize primary sequence space.
"""
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from tqdm import tqdm 

seq_path = "figures/figure3/primary_sequences.fasta"
names, seqs = [], []
with open(seq_path, "r") as fo:
    for line in fo:
        if line.startswith(">"): names.append(line.strip()[1:].split("|"))
        else: seqs.append(line.strip().split("|"))

def get_kmers(seqs, k):
    kmers = []
    for seq in seqs:
        for i in range(len(seq) - k + 1): kmers.append(" -> ".join(seq [ i : i + k ]))
    return kmers

# Mine for all unique kmers with k=3.
kmers3 = list(set(get_kmers(seqs, k=3)))
print(len(kmers3))

# Create dictionary of all unique kmers assigned to unique integer in a single dict.
kmers = {}
for i, kmer in enumerate(kmers3): kmers[kmer] = i

def kmer_fingerprint(seq):
    fp = np.zeros(len(kmers))
    fp_kmers = get_kmers([seq], k=3)
    for kmer in fp_kmers: fp[kmers[kmer]] += 1
    return fp

# Create fingerprint for each sequence.
fps = []
for seq in tqdm(seqs): fps.append(kmer_fingerprint(seq))
X = np.array(fps)
print("Fingerprints shape:", X.shape)

# Do PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X)
importances = pca.explained_variance_ratio_

# Plot PCA
plt.figure(figsize=(8, 8))
plt.scatter(X_pca[:, 0], X_pca[:, 1], c="black", s=10, alpha=0.9, label="Primary sequences (n=845)")
plt.xlabel("PC1 ({:.2f}%)".format(importances[0] * 100))
plt.ylabel("PC2 ({:.2f}%)".format(importances[1] * 100))
plt.xticks([])
plt.yticks([])
plt.gca().axes.xaxis.set_ticklabels([])
plt.gca().axes.yaxis.set_ticklabels([])
plt.gca().axes.xaxis.set_ticks([])
plt.gca().axes.yaxis.set_ticks([])
plt.legend() 
plt.tight_layout()
plt.savefig("figures/figure3/seqspace_pca.png", dpi=300)

# Per PC, get the 10 most important kmers.
loadings = pca.components_
for i in range(2):
    print("PC{}:".format(i + 1))
    for j in np.argsort(loadings[i])[::-1][:3]:
        print(kmers3[j], loadings[i][j])
