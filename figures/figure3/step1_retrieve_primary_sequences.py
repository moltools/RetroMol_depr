"""
Retrieve primary sequences.
"""
import argparse
import glob 
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm 
from rdkit import Chem 

from retromol.parsing import Result 

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, required=True, help="Path to dir with RetroMol JSON results.")
args = parser.parse_args()

# Define output path.
output_path = "figures/figure3/primary_sequences.fasta"

# Read RetroMol results JSONs.
filepaths = glob.glob(args.input + "/*.json")

# Parse results.
empty = 0
out_fo = open(output_path, "w")
for filepath in tqdm(filepaths, desc="Reading RetroMol results"):
    result = Result.from_json(filepath)
    if result.success and result.score == 0:
        seq = result.biosynthetic_seq
        if len(seq):
            seq = [x[1] for x in seq[0]]
            smi = Chem.MolToSmiles(result.substrate)
            out_fo.write(f">{result.name}|{smi}\n")
            out_fo.write("|".join(seq) + "\n")
        else:
            empty += 1
out_fo.close()
print(f"Number of empty sequences: {empty}")
