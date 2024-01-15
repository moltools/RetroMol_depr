#!/usr/bin/env python3
"""
Description:    Parse compounds from NPAtlas database JSON download.
Usage:          python3 parse_npatlas.py -i <path to NPAtlas JSON download> -o <path to output file>
"""
import argparse 
import json
from collections import Counter

def cli() -> argparse.Namespace:
    """
    Command line interface.
    
    :return:    Parsed command line arguments.
    :rtype:     argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True, type=str, help="path to NPAtlas JSON download.")
    parser.add_argument("-o", required=True, type=str, help="path to output file.")
    return parser.parse_args()

def main() -> None:
    """
    Main function.
    """
    args = cli()

    # Load JSON file.
    data = json.load(open(args.i))
    print("Number of compounds in dataset:", len(data))

    # Setup output file.
    out_fo = open(args.o, "w")
    out_fo.write("npaid\tsmiles\n")

    # Counters for classification types.
    pathway_classification_types = Counter()
    superclass_classification_types = Counter()
    class_classification_types = Counter()
    total_records = len(data) 
    parsed_records = 0

    # Iterate over each record in the JSON file.
    for record in data:
        try:
            npaid = record["npaid"]
            smiles = record["smiles"]
            pathway_classification = record["npclassifier"]["pathway_results"]
            superclass_classification = record["npclassifier"]["superclass_results"]
            class_classification = record["npclassifier"]["class_results"]

            for classification in pathway_classification:
                pathway_classification_types[classification] += 1

            # Filter on pathway classification.
            allowed_pathway_classifications = [
                "Polyketides", 
                "Amino acids and Peptides"
            ]

            # Filter on superclass classification.
            allowed_superclass_classifications = [
                "Oligopeptides", 
                "Macrolides", 
                "Cyclic polyketides", 
                "Linear polyketides", 
                "Small peptides", 
                "β-lactams", 
                "γ-lactam-β-lactones"
            ]

            # Filter on class classification.
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

            # Only write if all classifications are allowed.
            if (
                all([x in allowed_pathway_classifications for x in pathway_classification])
                and all([x in allowed_superclass_classifications for x in superclass_classification])
                and all([x in allowed_class_classifications for x in class_classification])
            ):

                for classification in superclass_classification: 
                    superclass_classification_types[classification] += 1

                for classification in class_classification: 
                    class_classification_types[classification] += 1

                out_fo.write(f"{npaid}\t{smiles}\n")
                parsed_records += 1

        except Exception:
            # Skip the record if it cannot be parsed.
            continue

        padding = len(str(total_records))
        print(f"{parsed_records}".zfill(padding) + "/" + f"{total_records}".zfill(padding), end="\r")

    out_fo.close()

    # Print most common pathway classification types.
    print("\nMost common pathway classification types in NPAtlas:")
    for i, (classification, count) in enumerate(pathway_classification_types.most_common()):
        print(f"{i}".zfill(2) + " " + f"{classification}: {count}")

    # Print most common superclass classification types.
    print("\nMost common superclass classification types in NPAtlas (after filtering on pathway classification):")
    for i, (classification, count) in enumerate(superclass_classification_types.most_common()):
        print(f"{i}".zfill(2) + " " + f"{classification}: {count}")

    # Print most common class classification types.
    print("\nMost common class classification types in NPAtlas (after filtering on pathway classification):")
    for i, (classification, count) in enumerate(class_classification_types.most_common()):
        print(f"{i}".zfill(2) + " " + f"{classification}: {count}")

    # Print number of parsed records.
    print(f"\nNumber of parsed records: {parsed_records}")