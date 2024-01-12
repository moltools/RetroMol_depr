"""
Parse polyketide, non-ribosomal peptide, RiPP, and hybrid natural products from NPAtlas. 
"""
import json
from collections import Counter
from tqdm import tqdm 

path_to_data = "npatlas/NPAtlas_download_20231129.json"
path_to_output = "npatlas/parsed_compounds.tsv"

data = json.load(open(path_to_data))
print("Number of compounds in dataset:", len(data)) # Expected: 33372 compounds

out_fo = open(path_to_output, "w")
out_fo.write("npaid\tsmiles\n")

pathway_classification_types = Counter()
superclass_classification_types = Counter()
class_classification_types = Counter()
parsed_records = 0

for record in tqdm(data, desc="Parsing records from NPAtlas"):
    try:
        npaid = record["npaid"]
        smiles = record["smiles"]
        pathway_classification = record["npclassifier"]["pathway_results"]
        superclass_classification = record["npclassifier"]["superclass_results"]
        class_classification = record["npclassifier"]["class_results"]

        for classification in pathway_classification:
            pathway_classification_types[classification] += 1

        allowed_pathway_classifications = [
            "Polyketides", 
            "Amino acids and Peptides"
        ]
        allowed_superclass_classifications = [
            "Oligopeptides", 
            "Macrolides", 
            "Cyclic polyketides", 
            "Linear polyketides", 
            "Small peptides", 
            "β-lactams", 
            "γ-lactam-β-lactones"
        ]
        allowed_class_classifications = [
            "Cyclic peptides", 
            "Depsipeptides", 
            "Linear peptides", 
            "Lipopeptides", 
            "Open-chain polyketides", 
            "Macrolide lactones", 
            "Polyene macrolides", 
            "Microcystins", 
            "Ansa macrolides", 
            "Dipeptides", 
            "Tripeptides", 
            "Macrolide lactams", 
            "Simple cyclic polyketides", 
            "Epothilones", 
            "Lactam bearing macrolide lactones", 
            "Salinosporamides"
        ]

        if (
            all([x in allowed_pathway_classifications for x in pathway_classification])
            and all([x in allowed_superclass_classifications for x in superclass_classification])
            and all([x in allowed_class_classifications for x in class_classification])
        ):

            for classification in superclass_classification: superclass_classification_types[classification] += 1
            for classification in class_classification: class_classification_types[classification] += 1

            out_fo.write(f"{npaid}\t{smiles}\n")
            parsed_records += 1

    except Exception:
        continue

out_fo.close()

print("\nMost common pathway classification types in NPAtlas:")
for i, (classification, count) in enumerate(pathway_classification_types.most_common()):
    print(f"{i}".zfill(2) + " " + f"{classification}: {count}")

print("\nMost common superclass classification types in NPAtlas (after filtering on pathway classification):")
for i, (classification, count) in enumerate(superclass_classification_types.most_common()):
    print(f"{i}".zfill(2) + " " + f"{classification}: {count}")

print("\nMost common class classification types in NPAtlas (after filtering on pathway classification):")
for i, (classification, count) in enumerate(class_classification_types.most_common()):
    print(f"{i}".zfill(2) + " " + f"{classification}: {count}")

print(f"\nNumber of parsed records: {parsed_records}") # Expected: 16603 compounds
