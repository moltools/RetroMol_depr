#!/usr/bin/env python3
"""This module contains the command line interface for the RetroMol package."""
import argparse
import json
import logging
import multiprocessing as mp
import os
import typing as ty

from rdkit import RDLogger

from retromol.chem import Molecule, MolecularPattern, ReactionRule

from retromol.helpers import timeout
from retromol.parsing import (
    Result,
    parse_reaction_rules,
    parse_molecular_patterns,
    parse_mol,
)

# RDLogger.DisableLog("rdApp.*")

logger = logging.getLogger(__name__)


def cli() -> argparse.Namespace:
    """Parse command line arguments.

    :return: The parsed command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit.",
    )
    parser.add_argument(
        "-r",
        "--reactions",
        type=str,
        required=False,
        default=os.path.join(os.path.dirname(os.path.dirname(__file__)), "tests", "fixtures", "reactions.json"),
        help="Path to JSON file containing reaction.",
    )
    parser.add_argument(
        "-m",
        "--monomers",
        type=str,
        required=False,
        default=os.path.join(os.path.dirname(os.path.dirname(__file__)), "tests", "fixtures", "monomers.json"),
        help="Path to JSON file containing monomer patterns.",
    )
    parser.add_argument(
        "-l",
        "--logger-level",
        type=str,
        default="INFO",
        help="Logger verbosity level.",
    )

    subparsers = parser.add_subparsers(dest="mode", required=True)
    subparser_single_mode = subparsers.add_parser("single")
    subparser_batch_mode = subparsers.add_parser("batch")

    subparser_single_mode.add_argument(
        "-i", "--input", type=str, required=True, help="Input SMILES (between quotes)."
    )
    subparser_single_mode.add_argument(
        "-o", "--output", type=str, required=True, help="Path to new JSON output file."
    )

    subparser_batch_mode.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Path to input CSV or TSV file as 'name,SMILES' per line.",
    )
    subparser_batch_mode.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Path to new directory for JSON output files.",
    )
    subparser_batch_mode.add_argument(
        "-n",
        "--nproc",
        type=int,
        default=mp.cpu_count(),
        help="Number of processes to use.",
    )
    subparser_batch_mode.add_argument(
        "--header", action="store_true", help="Input file contains a header line."
    )

    return parser.parse_args()


Record = ty.Tuple[Molecule, ty.List[ReactionRule], ty.List[MolecularPattern]]


@timeout(10)
def parse_mol_timed(record: Record) -> Result:
    """Parse a molecule with a timeout.

    :param record: The record containing the molecule, reaction rules, and
        molecular patterns.
    :type record: Record
    :return: The result of the parsing.
    :rtype: Result
    """
    logger = logging.getLogger(__name__)
    try:
        mol, reactions, monomers = record
        return parse_mol(mol, reactions, monomers)
    except Exception as e:
        logger.error(f"Failed to parse {mol.name} and raised {e.__class__.__name__}: {e}")
        return Result(mol.name, mol.compiled, False)


def main() -> None:
    """
    Driver function.
    """
    # Disable RDKit logging.
    RDLogger.DisableLog("rdApp.*")

    # Parse command line arguments.
    args = cli()

    # Set logger verbosity.
    logging.basicConfig(level=args.logger_level)

    # Parse reaction rules and molecular patterns.
    reactions_src = json.load(open(args.reactions, "r", encoding="utf-8"))
    reactions = parse_reaction_rules(reactions_src)

    monomers_src = json.load(open(args.monomers, "r", encoding="utf-8"))
    monomers = parse_molecular_patterns(monomers_src)
    
    logger.debug(f"Parsed {len(reactions)} reactions and {len(monomers)} monomers.")

    # Parse molecule in single mode.
    if args.mode == "single":
        mol = Molecule("input", args.input)
        record = (mol, reactions, monomers)
        result = parse_mol_timed(record)
        with open(args.output, "w", encoding="utf-8") as fo:
            fo.write(result.to_json())
        logger.debug(f"Processed single molecule.")
        exit(0)

    # Parse a batch of molecules.
    elif args.mode == "batch":

        # Parse input file.
        records = []
        with open(args.input, "r", encoding="utf-8") as fo:
            # Check file extension.
            sep = {"csv": ",", "tsv": "\t"}.get(args.input.split(".")[-1])
            if args.header:
                next(fo)
            for i, line in enumerate(fo):
                name, smiles = line.strip().split(sep)
                mol = Molecule(name, smiles)
                records.append((mol, reactions, monomers))

                # Report on progress.
                print(f"Added record {i+1} to queue...", end="\r")

        # Parse molecules in parallel.
        nproc = min(args.nproc, mp.cpu_count())
        print(f"\nUsing {nproc} processes...")
        with mp.Pool(processes=nproc) as pool:
            for i, result in enumerate(pool.imap_unordered(parse_mol_timed, records)):
                out_path = os.path.join(args.output, f"{result.identifier}.json")
                with open(out_path, "w", encoding="utf-8") as fo:
                    fo.write(result.to_json())

                # Report on progress.
                percentage_done = round((i + 1) / len(records) * 100, 2)
                print(f"Processed {percentage_done}% of records...", end="\r")
        logger.info(f"Processed {len(records)} records.")
        exit(0)

    # Invalid mode.
    else:
        print(f"Invalid mode: {args.mode}")
        exit(1)


if __name__ == "__main__":
    main()
