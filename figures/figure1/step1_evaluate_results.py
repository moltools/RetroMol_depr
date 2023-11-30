"""
Evaluate RetroMol results with a Sankey graph and scores histogram.

Dependencies:
    - matplotlib
    - numpy
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

# Read RetroMol results JSONs.
filepaths = glob.glob(args.input + "/*.json")

# Parse results.
scores = []
scores_as_percentages = []
number_of_atoms = []
total, failed, success = 0, 0, 0

for filepath in tqdm(filepaths, desc="Reading RetroMol results"):
    result = Result.from_json(filepath)
    
    total += 1
    if result.success:
        success += 1

        # Get score.
        scores.append(result.score)

        # Calculate score as percentage.
        total_atoms = result.substrate.GetNumAtoms()
        number_of_atoms.append(total_atoms)

        missing_atoms = result.score
        score_as_percentage = missing_atoms / total_atoms * 100
        score_as_percentage = round(score_as_percentage)
        scores_as_percentages.append(score_as_percentage)

    else:
        failed += 1

# Print. 
print(f"\nSuccess: {success}")
print(f"Failed: {failed}")
print(f"Total: {total}")

# Print number of compounds with score 0.
print(f"\nNumber of compounds with score 0: {scores.count(0)}")

# Print number of compounds with unaccounted score of 100:
print(f"Number of compounds with score 100: {scores_as_percentages.count(100)}")

# Get max score and make histogram for every int between 0 and that max score.
max_score = max(scores)

# Make histogram for scores.
hist, bins = np.histogram(scores, bins=max_score)
center = (bins[:-1] + bins[1:]) / 2
fig, ax = plt.subplots()
ax.bar(center, hist, align="center")
ax.set_xlabel("Score (num. of atoms unaccounted for)")
ax.set_ylabel("Frequency (num. of compounds)")
ax.set_title("Monomer coverage")
plt.savefig("figures/figure1/histogram_scores_as_counts.png", dpi=300)

# Make histogram for scores as percentages.
hist, bins = np.histogram(scores_as_percentages, bins=100)
center = (bins[:-1] + bins[1:]) / 2
fig, ax = plt.subplots()
ax.bar(center, hist, align="center")
ax.set_xlabel("Score (perc. of atoms unaccounted for)")
ax.set_ylabel("Frequency (num. of compounds)")
ax.set_title("Monomer coverage")
ax.xaxis.set_major_formatter('{x:.0f}%')
plt.savefig("figures/figure1/histogram_scores_as_percentages.png", dpi=300)

# Make histogram for number of atoms.
hist, bins = np.histogram(number_of_atoms, bins=100)
center = (bins[:-1] + bins[1:]) / 2
fig, ax = plt.subplots()
ax.bar(center, hist, align="center")
ax.set_xlabel("Number of atoms")
ax.set_ylabel("Frequency (num. of compounds)")
ax.set_title("Size of compounds")
plt.savefig("figures/figure1/histogram_number_of_atoms.png", dpi=300)
