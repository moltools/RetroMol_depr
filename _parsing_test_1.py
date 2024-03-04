import argparse 
import typing as ty
import json

from retromol.chem import Molecule, MolecularPattern, ReactionRule
from retromol.helpers import TimeoutError, timeout
from retromol.parsing import Result, parse_reaction_rules, parse_molecular_patterns, parse_mol 
from retromol_sequencing.alignment import ModuleSequence, parse_primary_sequence, MultipleSequenceAlignment
from retromol_sequencing.sequencing import parse_modular_natural_product

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-r", "--rules", type=str, required=True, help="Path to JSON file containing rules.")   
    return parser.parse_args()

def main() -> None:
    args = cli()

    rules = json.load(open(args.rules, "r"))
    reactions = parse_reaction_rules(json.dumps(rules["reactions"]))
    monomers = parse_molecular_patterns(json.dumps(rules["monomers"]))

    # erythromycin
    smiles = r"CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O"

    # anguinomycin
    # smiles = r"CC/C(=C/C(C)C/C=C/C(=C/C(C)C(=O)C(C)C(C(C)C/C(=C/C)/C)O)/C)/C=C/C1CC=CC(=O)O1"

    mol = Molecule("molecule", smiles)

    result = parse_mol(mol, reactions, monomers)
    sequences = parse_modular_natural_product(result)

    for i, s in enumerate(sequences):
        for motif in s:
            print(i, motif)

if __name__ == "__main__":
    main()