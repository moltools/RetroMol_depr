#!/usr/bin/env python3
"""
Description:    RetroMol command line interface.
Usage:          python3 main.py -r <path/to/rules.json> <mode> <args>
"""
import argparse
import json 
import multiprocessing as mp
import os 
import typing as ty 

from rdkit import RDLogger
from tqdm import tqdm 

from retromol.chem import Molecule, MolecularPattern, ReactionRule
from retromol.helpers import TimeoutError, timeout
from retromol.parsing import Result, parse_reaction_rules, parse_molecular_patterns, parse_mol 

def cli() -> argparse.Namespace:
    """
    Command line interface.

    :returns: Command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help="Show this help message and exit.")
    parser.add_argument("-r", "--rules", type=str, required=True, help="Path to JSON file containing rules.")

    subparsers = parser.add_subparsers(dest="mode", required=True)
    subparser_single_mode = subparsers.add_parser("single")
    subparser_batch_mode = subparsers.add_parser("batch")

    subparser_single_mode.add_argument("-i", "--input", type=str, required=True, help="Input SMILES (between quotes).")
    subparser_single_mode.add_argument("-o", "--output", type=str, required=True, help="Path to new JSON output file.")

    subparser_batch_mode.add_argument("-i", "--input", type=str, required=True, help="Path to input TSV file as 'name\tSMILES'.")
    subparser_batch_mode.add_argument("-o", "--output", type=str, required=True, help="Path to new directory for JSON output files.")
    subparser_batch_mode.add_argument("-n", "--nproc", type=int, default=mp.cpu_count(), help="Number of processes to use.")
    subparser_batch_mode.add_argument("--header", action="store_true", help="Input file contains a header line.")

    return parser.parse_args()

Record = ty.Tuple[Molecule, ty.List[ReactionRule], ty.List[MolecularPattern]]

@timeout(10)
def parse_mol_timed(record: Record) -> Result:
    """
    Parse molecule.
    
    :param Record record: Record containing molecule, reaction rules, and molecular patterns.
    :returns: Result object.
    :rtype: Result
    """
    try:
        mol, reactions, monomers = record
        return parse_mol(mol, reactions, monomers)
    except Exception:
        return Result(mol.name, mol.compiled, False)

def main() -> None:
    """
    Driver function.
    """
    # Disable RDKit logging.
    RDLogger.DisableLog("rdApp.*")

    # Parse command line arguments.
    args = cli()

    # Parse reaction rules and molecular patterns.
    rules = json.load(open(args.rules, "r"))
    reactions = parse_reaction_rules(json.dumps(rules["reactions"]))
    monomers = parse_molecular_patterns(json.dumps(rules["monomers"]))

    # Parse molecule in single mode.
    if args.mode == "single":
        mol = Molecule("input", args.input)
        record = (mol, reactions, monomers)
        result = parse_mol_timed(record)
        with open(args.output, "w") as fo: fo.write(result.to_json())

    # Parse a batch of molecules.
    elif args.mode == "batch":
        
        # Parse input file.
        records = []
        with open(args.input, "r") as fo:
            if args.header: next(fo)
            for line in tqdm(fo, desc="Parsing records from input file"):
                name, smiles = line.strip().split("\t")
                mol = Molecule(name, smiles)
                records.append((mol, reactions, monomers))

        # Parse molecules in parallel.
        nproc = min(args.nproc, mp.cpu_count())
        with mp.Pool(processes=nproc) as pool:
            for result in tqdm(pool.imap_unordered(parse_mol_timed, records)):
                out_path = os.path.join(args.output, f"{result.name}.json")
                with open(out_path, "w") as fo: fo.write(result.to_json())

    # Invalid mode.
    else:
        print(f"Invalid mode: {args.mode}")

    exit("Done.")

if __name__ == "__main__":
    main()