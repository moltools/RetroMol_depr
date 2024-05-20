#!/usr/bin/env python3
"""Annotate biosynthetic sequences in the Neo4j database."""
import argparse
import logging
import os
import re
import typing as ty
import json

from neo4j import GraphDatabase
from rdkit import RDLogger
from tqdm import tqdm

from retromol.retrosynthesis.parsing import Result
from retromol.retrosynthesis.sequencing import parse_modular_natural_product

RDLogger.DisableLog("rdApp.*")


def cli() -> argparse.Namespace:
    """Parse command-line arguments.

    :return: Command-line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--results", type=str, required=True, help="Path to dir containing JSON results."
    )
    parser.add_argument("--port", default=7687, type=int, help="Neo4j port.")
    parser.add_argument(
        "--authentication",
        default=None,
        nargs=2,
        help="Neo4j authentication as '<username> <password>'.",
    )
    return parser.parse_args()


def parse_motif(i: int, motif: str) -> str:
    """Parse a motif into a Cypher query.

    :param i: Index of motif.
    :type i: int
    :param motif: Motif to parse.
    :type motif: str
    :return: Cypher query.
    :rtype: str
    """
    if match := re.match(r"polyketide\|([A-D])(\d{1,2})", motif):
        accessory_domains = match.group(1)
        decoration_type = int(match.group(2))

        return (
            f"(u{i + 1}:Motif {{"
            f"identifier: 'polyketide', "
            f"accessory_domains: '{accessory_domains}', "
            f"decoration_type: '{decoration_type}'"
            f"}})"
        )

    elif match := re.match(r"peptide\|(\w+)\|(.+)", motif):
        source = match.group(1)
        cid = match.group(2)

        return (
            f"(u{i + 1}:Motif {{"
            f"identifier: 'peptide', "
            f"source: '{source}', "
            f"cid: '{cid}'"
            f"}})"
        )

    else:
        raise ValueError(f"Unknown motif: {motif}")


def main() -> None:
    """Driver function."""
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    args = cli()

    # Get all paths to JSON files in results directory.
    paths = [
        os.path.join(args.results, path)
        for path in os.listdir(args.results)
        if path.endswith(".json")
    ]
    logger.info(f"Found {len(paths)} JSON files in {args.results}.")

    if args.authentication:
        db = GraphDatabase.driver(f"bolt://localhost:{args.port}", auth=tuple(args.authentication))
    else:
        db = GraphDatabase.driver(f"bolt://localhost:{args.port}")

    with db.session() as session:

        # Iterate over all JSON files.
        for path in tqdm(paths):

            result = Result.from_json(path)

            if result.success is False:
                continue

            identifier = result.identifier  # We assume that the identifier is the NPAtlas ID.
            applied_reactions = result.applied_reactions
            sequences = result.sequences

            # Add biosynthetic fingerprint to graph.
            for accession, sequence in enumerate(sequences):

                # Create biosynthetic fingerprint node sequence.
                query_begin = """
                CREATE (s:PrimarySequence {identifier: $identifier, accession: $accession, applied_reactions: $applied_reactions, sequences: $sequences})-[:START]->
                """
                query_end = "-[:NEXT]->".join(
                    [parse_motif(i, motif) for i, motif in enumerate(sequence["motif_code"])]
                )
                session.run(
                    query_begin + query_end,
                    identifier=identifier,
                    accession=accession,
                    applied_reactions=applied_reactions,
                    sequences=json.dumps(sequences)
                )

                # Add link between biosynthetic fingerprint and compound.
                query = """
                MATCH (c:Compound {identifier: $identifier})
                MATCH (s:PrimarySequence {identifier: $identifier})
                CREATE (c)-[:HAS_PRIMARY_SEQUENCE]->(s)
                """
                session.run(query, identifier=identifier)

    db.close()


if __name__ == "__main__":
    main()
