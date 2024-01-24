#!/usr/bin/env python3
"""
Description:    This script retrieves the primary sequence from a monomer graph.
Usage:          python3 retrieve_primary_sequence.py -i <path/to/retromol/output/json>
"""
import argparse
import glob
import joblib
import json 
import os
import re 
import typing as ty

from rdkit import Chem 

from retromol.parsing import Result, parse_molecular_patterns
from retromol_sequencing.fingerprint import get_amino_acid_fingerprint, amino_acid_class_to_label
from retromol_sequencing.primary_sequence import resolve_biosynthetic_sequence

def cli() -> argparse.Namespace:
    """
    Command line interface.

    :returns: Command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="Path to dir with result JSON files.")
    parser.add_argument("-m", "--model", type=str, required=True, help="Path to model file for classifying monomers.")
    parser.add_argument("-r", "--rules", type=str, required=True, help="Input rules JSON file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to output .seq file.")
    return parser.parse_args()

def is_polyketide(name: str, motif: Chem.Mol) -> bool:
    """
    Check if motif is a polyketide.

    :param str name: Name of motif.
    :param Chem.Mol motif: Motif.
    :returns: True if motif is a polyketide.
    :rtype: bool
    """
    if re.match(r"^[A-D]\d+$", name):
        return True
    
    return False

def is_amino_acid(name: str, motif: Chem.Mol) -> bool:
    """
    Check if motif is an amino acid.

    :param str name: Name of motif.
    :param Chem.Mol motif: Motif.
    :returns: True if motif is an amino acid.
    :rtype: bool
    """
    smarts_strings = [
        r"[NH1,NH2][C][C](=[O])[O]", # Alpha amino acid
        r"[NH1,NH2][C][C][C](=[O])[O]" # Beta amino acid
    ]

    if any([motif.HasSubstructMatch(Chem.MolFromSmarts(x)) for x in smarts_strings]):
        return True
    
    return False

def main() -> None:
    """
    Main function.
    """
    # Parse command line arguments.
    args = cli()

    # Get paths to all JSON files in input directory.
    paths = glob.glob(os.path.join(args.input, "*.json"))

    # Load reaction rules.
    rules = json.load(open(args.rules, "r"))
    monomers = parse_molecular_patterns(json.dumps(rules["monomers"]))

    # Load model.
    model = joblib.load(args.model)

    # Define function for checking if a motif is a primary motif.
    def is_primary_motif(motif: Chem.Mol) -> ty.Optional[str]:
            """
            Check if motif is a primary motif.

            :param Chem.Mol motif: Motif.
            :returns: Primary motif identity or None.
            :rtype: ty.Optional[str]
            """
            for monomer in monomers:
                if motif.HasSubstructMatch(monomer.compiled):
                    name = monomer.name

                    if is_polyketide(name, motif):
                        return name 
                    elif is_amino_acid(name, motif):
                        return name
                
            return False
    
    # Classify monomer.
    def classify(name: str, mol: Chem.Mol) -> str:
        """
        Classify monomer.
        
        :param str name: Name of monomer.
        :param Chem.Mol mol: Molecule representing monomer.
        :returns: Classification.
        :rtype: str
        """
        if is_polyketide(name, mol):
            return name[0] # Any of A, B, C, D.
        else:
            fp = get_amino_acid_fingerprint(mol)
            pred = model.predict([fp])[0]
            label = amino_acid_class_to_label(pred)
            return label

    # Iterate over all JSON files.
    out_file = open(args.output, "w")

    for i, path in enumerate(paths):
        # Load result.
        result = Result.from_json(path)

        if result.success == True:
            # Resolve primary sequence.
            seq = resolve_biosynthetic_sequence(result, is_primary_motif)

            # Process primary sequence further if it is long enough.
            if len(seq) >= 2:

                # Compile header.
                name = result.name 
                score = result.score
                smiles = Chem.MolToSmiles(result.substrate)
                header = f">{name}|{score}|{smiles}\n"
                out_file.write(header)

                # Get motif classification for each monomer based on chemical structure.
                seq = [(t[0], classify(t[0], t[1])) for t in seq]
                
                # Split tuples into two sequences.
                seq1, seq2 = zip(*seq)
                out_file.write("|".join(seq1) + "\n")
                out_file.write("|".join(seq2) + "\n")
                out_file.flush()

        # Print progress.
        perc = round((i + 1) / len(paths) * 100, 2)
        print(f"{perc}%", end="\r")

    # Write results to JSON file.
    out_file.close()

    exit(0)

if __name__ == "__main__":
    main()