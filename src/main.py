#!/usr/bin/env python3
import argparse
import errno 
import functools
import json 
import multiprocessing as mp
import os 
import signal
import typing as ty 

from rdkit import RDLogger
from tqdm import tqdm 

from retromol.chem import Molecule, MolecularPattern, ReactionRule
from retromol.parsing import Result, parse_reaction_rules, parse_molecular_patterns, parse_mol 

class TimeoutError(Exception):
    pass 

def timeout(seconds: int = 5, error_message: str = os.strerror(errno.ETIME)) -> ty.Callable:
    """
    Decorator to raise a TimeoutError when runtime exceeds the specified time.

    Parameters
    ----------
    seconds : int, optional
        Timeout in seconds, by default 5.
    error_message : str, optional
        Error message to be raised, by default os.strerror(errno.ETIME).
    
    Note: Timer only works on UNIX systems; also not thread-safe.
    """
    def decorator(func: ty.Callable) -> ty.Callable:

        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)
        
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try: 
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result
        
        return wrapper
    
    return decorator

def cli() -> argparse.Namespace:
    """
    Command line interface.
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

@timeout(5)
def parse_mol_timed(record: Record) -> Result:
    """
    Parse molecule.
    
    Parameters
    ----------
    record : Record
        Tuple containing molecule, reaction rules, and motif units.
    
    Returns
    -------
    result : Result
        Result object.
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
    RDLogger.DisableLog("rdApp.*")
    args = cli()

    rules = json.load(open(args.rules, "r"))
    reactions = parse_reaction_rules(json.dumps(rules["reactions"]))
    monomers = parse_molecular_patterns(json.dumps(rules["monomers"]))

    if args.mode == "single":
        mol = Molecule("input", args.input)
        record = (mol, reactions, monomers)
        result = parse_mol_timed(record)
        with open(args.output, "w") as fo: fo.write(result.to_json())

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

    else:
        print(f"Invalid mode: {args.mode}")

    exit("Done.")

if __name__ == "__main__":
    main()