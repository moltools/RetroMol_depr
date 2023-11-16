#!/usr/bin/env python3
import argparse 
import os

from rdkit import RDLogger, Chem
from rdkit.Chem import AllChem
RDLogger.DisableLog('rdApp.*')

from retromol.chem import draw_molecule
from retromol.parsing import Result, generate_structures, parse_reaction_rules, parse_molecular_patterns

def cli() -> argparse.Namespace:   
    """
    Command-line interface.
    """
    parser = argparse.ArgumentParser()
    # parser.add_argument("-i", "--input", type=str, required=True, help="Path to JSON file containing monomer graph.")
    parser.add_argument("-rr", "--reaction-rules", type=str, required=True, help="Path to reaction rules JSON file.")
    parser.add_argument("-cu", "--core-units", type=str, required=True, help="Path to core units JSON file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to output dir.")
    return parser.parse_args()

def main() -> None:
    """
    Driver code.
    """
    args = cli()
    reactions = parse_reaction_rules(args.reaction_rules)
    core_units = parse_molecular_patterns(args.core_units, as_smiles=True)
    # result = Result.from_json(args.input)
    # recipe = result.get_biosynthetic_recipe()
    # mols = generate_structures(recipe, reactions, core_units)
    # print(mols)

    reactions = {r.name: r.forward_compiled for r in reactions}
    pks = reactions["PKS"]

    def selected(name: str) -> bool:
        return name[0] in ["A", "B", "C", "D"] and name[1:].isdigit()

    core_units = {c.name: c.compiled for c in core_units if selected(c.name)}
    
    seq = ["B7", "B2", "A2", "D7", "B2", "B2"]
    curr = core_units[seq[0]]
    for x in seq[1:]:
        x = core_units[x]
        reactants = [curr, x]
        results = pks.RunReactants(reactants)
        if len(results) == 1:
            curr = results[0][0]
            Chem.SanitizeMol(curr)
        else:
            print("Failed")
            break

    smirks = r"[OH][S][C:3]>>[C][C][C:3]"
    rxn = AllChem.ReactionFromSmarts(smirks)
    results = rxn.RunReactants([curr])
    product = results[0][0]
    Chem.SanitizeMol(product)
    curr = product

    cycl = reactions["macrocyclization"]
    results = cycl.RunReactants([curr])
    for i, result in enumerate(results):
        product = result[0]
        Chem.SanitizeMol(product)
        
        smiles = Chem.MolToSmiles(product)
        print(smiles)

        out_path = os.path.join(args.output, f"option_{i}.png")
        draw_molecule(product, out_path)

    exit(0)

if __name__ == "__main__":
    main()