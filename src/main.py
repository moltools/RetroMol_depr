#!/usr/bin/env python3
import argparse
import errno 
import functools
import os 
import signal
import typing as ty 

from retromol.chem import Molecule, MolecularPattern, ReactionRule
from retromol.parsing import (
    Result, 
    parse_reaction_rules, 
    parse_molecular_patterns, 
    parse_mol 
)

class TimeoutError(Exception):
    pass 

def timeout(
    seconds: int = 5, 
    error_message: str = os.strerror(errno.ETIME)
) -> ty.Callable:
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
    parser.add_argument("-rr", "--reaction-rules", type=str, required=True, help="Path to file containing reaction rules.")
    parser.add_argument("-cu", "--core-units", type=str, required=True, help="Path to file containing core motif units.")
    parser.add_argument("-ou", "--other-units", type=str, required=True, help="Path to file containing other monomer units.")

    subparsers = parser.add_subparsers(dest="mode", required=True)
    subparser_single_mode = subparsers.add_parser("single")
    subparser_batch_mode = subparsers.add_parser("batch")

    subparser_single_mode.add_argument("-i", "--input", type=str, required=True, help="Input SMILES.")

    return parser.parse_args()

@timeout()
def parse_mol_timed(
    mol: Molecule, 
    reactions: ty.List[ReactionRule],
    core_units: ty.List[MolecularPattern],
    other_units: ty.List[MolecularPattern]
) -> Result:
    """
    Parse molecule.

    Parameters
    ----------
    mol : Molecule
        Molecule.
    reactions : ty.List[ReactionRule]
        List of reaction rules.
    core_units : ty.List[MolecularPattern]
        List of motif units.
    other_units : ty.List[MolecularPattern]
        List of other monomer units.
    
    Returns
    -------
    result : Result
        Result object.
    """
    return parse_mol(mol, reactions, core_units, other_units)

def main() -> None:
    """
    Driver function.
    """
    args = cli()

    rxn = parse_reaction_rules(args.reaction_rules)
    cus = parse_molecular_patterns(args.core_units)
    ous = parse_molecular_patterns(args.other_units)

    match args.mode:
        case "single":
            mol = Molecule("input", args.input)
            result = parse_mol_timed(mol, rxn, cus, ous)
            print(result.to_json())
        
        case "batch":
            print("Batch mode not implemented yet.")

        case _:
            print(f"Invalid mode: {args.mode}")

    exit(0)

if __name__ == "__main__":
    main()