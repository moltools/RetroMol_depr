#!/usr/bin/env python3
import argparse 
import random
from copy import deepcopy

from rdkit import RDLogger, Chem
from rdkit.Chem import AllChem
RDLogger.DisableLog('rdApp.*')

from retromol.parsing import parse_reaction_rules, parse_molecular_patterns

def cli() -> argparse.Namespace:   
    """
    Command-line interface.
    """
    parser = argparse.ArgumentParser()
    # parser.add_argument("-i", "--input", type=str, required=True, help="Path to JSON file containing monomer graph.")
    parser.add_argument("-rr", "--reaction-rules", type=str, required=True, help="Path to reaction rules JSON file.")
    parser.add_argument("-cu", "--core-units", type=str, required=True, help="Path to core units JSON file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to output file.")
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

    # Generate 10.000 unique sequences randomly of length n from core_units keys, with replacement:
    n = 1000
    size = 6
    seqs = set()
    while len(seqs) < n:
        seq = "|".join(random.choices(list(core_units.keys()), k=size))
        seqs.add(seq)

    mols = []
    for seq in seqs:
        try:
            seq = seq.split("|")
            curr = core_units[seq[0]]
            for x_name in seq[1:]:
                x = core_units[x_name]
                reactants = [curr, x]
                results = pks.RunReactants(reactants)
                if len(results) == 1:
                    curr = results[0][0]
                    Chem.SanitizeMol(curr)
                else:
                    raise(Exception("Reaction failed."))
        except Exception as e:
            print(e)
            continue

        mols.append((curr, seq))

    duplicate_check = set()
    cyclized_mols = []
    for mol in mols:
        curr = mol[0]
        cycl = reactions["macrocyclization"]
        results = cycl.RunReactants([curr])
        for i, result in enumerate(results):
            product = result[0]
            Chem.SanitizeMol(product)
            smi = Chem.MolToSmiles(product)
            if smi not in duplicate_check:
                duplicate_check.add(smi)
                cyclized_mols.append((product, mol[1]))
    
    mols += cyclized_mols

    smirks = r"[OH][S][C:3]>>[C][C][C:3]"
    rxn = AllChem.ReactionFromSmarts(smirks)
    capped_mols = []
    for mol, seq in mols:
        results = rxn.RunReactants([mol])
        product = results[0][0]
        Chem.SanitizeMol(product)
        capped_mols.append((product, seq))

    mols = capped_mols
    smis = [(Chem.MolToSmiles(mol), seq) for mol, seq in mols]

    with open (args.output, "w") as f:
        for smi, seq in smis:
            seq = "|".join(seq)
            f.write(f"{smi},{seq}\n")
    exit(0)

if __name__ == "__main__":
    main()