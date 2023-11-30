"""
Visualize frequent kmers.
"""
import glob 
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
from tqdm import tqdm 

from retromol.parsing import Result 

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

# Analayze k=2
kmers = Counter(get_kmers(seqs, k=2))
plt.figure(figsize=(8, 15))
plt.barh(range(20), [x[1] for x in kmers.most_common(20)])
plt.yticks(range(20), [x[0] for x in kmers.most_common(20)])
plt.xlabel("Frequency")
plt.tight_layout()
plt.savefig("figures/figure3/barplot_k2.png", dpi=300)

# Analayze k=3
kmers = Counter(get_kmers(seqs, k=3))
plt.figure(figsize=(8, 15))
plt.barh(range(20), [x[1] for x in kmers.most_common(20)])
plt.yticks(range(20), [x[0] for x in kmers.most_common(20)])
plt.xlabel("Frequency")
plt.tight_layout()
plt.savefig("figures/figure3/barplot_k3.png", dpi=300)

# Analayze k=4
kmers = Counter(get_kmers(seqs, k=4))
plt.figure(figsize=(8, 15))
plt.barh(range(20), [x[1] for x in kmers.most_common(20)])
plt.yticks(range(20), [x[0] for x in kmers.most_common(20)])
plt.xlabel("Frequency")
plt.tight_layout()
plt.savefig("figures/figure3/barplot_k4.png", dpi=300)

# Analayze k=5
kmers = Counter(get_kmers(seqs, k=5))
plt.figure(figsize=(8, 15))
plt.barh(range(20), [x[1] for x in kmers.most_common(20)])
plt.yticks(range(20), [x[0] for x in kmers.most_common(20)])
plt.xlabel("Frequency")
plt.tight_layout()
plt.savefig("figures/figure3/barplot_k5.png", dpi=300)

