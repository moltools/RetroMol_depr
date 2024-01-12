#!/usr/bin/env python3
import argparse 
import glob
import os 
from tqdm import tqdm
from rdkit import Chem 
from retromol.parsing import Result 
from retromol.drawing import draw_molecule

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Draw failed compounds")
    parser.add_argument("--input", type=str, required=True, help="Path to RetroMol results dir.")
    parser.add_argument("--output", type=str, required=True, help="Path to output dir.")
    return parser.parse_args()

def main() -> None:
    args = cli()

    # Read all json failes from input dir.
    files = glob.glob(args.input + "/*.json")
    print("Found {} json files.".format(len(files)))

    # Parse results.
    num_parsed = 0
    for file in tqdm(files):
        result = Result.from_json(file)
        if result.success and result.score > 0:
            num_parsed += 1

            npaid = result.name

            # Get SMILES.
            substrate = result.substrate
            path_out = os.path.join(args.output, "{}.png".format(npaid))
            draw_molecule(substrate, path_out)

    print("Parsed {} compounds.".format(num_parsed))

    exit(0)
    
if __name__ == "__main__":
    main()    

