"""
Parse polyketide, non-ribosomal peptide, RiPP, and hybrid natural products from NPAtlas. 
"""
import json
from collections import Counter
from tqdm import tqdm 

path_to_data = "npatlas/NPAtlas_download_20231129.json"
path_to_output = "npatlas/parsed_compounds_20231129.tsv"

data = json.load(open(path_to_data))
print("Number of compounds in dataset:", len(data)) # Expected: 33372 compounds

out_fo = open(path_to_output, "w")
out_fo.write("npaid\tsmiles\n")

pathway_classification_types = Counter()
parsed_records = 0

for record in tqdm(data, desc="Parsing records from NPAtlas"):
    try:
        npaid = record["npaid"]
        smiles = record["smiles"]
        pathway_classification = record["npclassifier"]["pathway_results"]

        for classification in pathway_classification:
            pathway_classification_types[classification] += 1

        allowed_pathway_classifications = ["Polyketides", "Amino acids and Peptides"]

        if all([x in allowed_pathway_classifications for x in pathway_classification]):
            out_fo.write(f"{npaid}\t{smiles}\n")
            parsed_records += 1

    except Exception:
        continue

out_fo.close()

print("\nMost common pathway classification types in NPAtlas:")
for i, (classification, count) in enumerate(pathway_classification_types.most_common()):
    print(f"{i}".zfill(2) + " " + f"{classification}: {count}")

print(f"\nNumber of parsed records: {parsed_records}") # Expected: 16603 compounds
