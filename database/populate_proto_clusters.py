#!/usr/bin/env python3
"""This script populates the Neo4j database with protoclusters."""
import argparse
import glob
import json
import logging
import re

from neo4j import GraphDatabase
from tqdm import tqdm

def cli() -> argparse.Namespace:
    """Parse command-line arguments.
    
    :return: Command-line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to directory containing parsed data.")
    parser.add_argument("--port", default=7687, type=int, help="Neo4j port.")
    parser.add_argument(
        "--authentication",
        default=None,
        nargs=2,
        help="Neo4j authentication as '<username> <password>'."
    )
    return parser.parse_args()


def parse_json(file_path: str):
    """
    """
    logger = logging.getLogger(__name__)

    with open(file_path, "r") as file:
        data = json.load(file)

    for data_item in data:
        name = data_item["name"]
        description = data_item["description"]
        protoclusters = data_item["proto_clusters"]

        logger.debug(f"Found {len(protoclusters)} protoclusters in {name}.")

        for protocluster in protoclusters:
            for protocluster_item in protocluster:
                
                sequence = []

                accession = protocluster_item["accession"]
                modules = protocluster_item["modules"]

                for module in modules:
                    hit_ids = set()
                    for component in module["components"]:
                        hit_id = component["domain"]["hit_id"]
                        hit_ids.add(hit_id)

                    if any([hit_id.startswith("PKS") for hit_id in hit_ids]):
                        has_kr = any([hit_id.startswith("PKS_KR") for hit_id in hit_ids])
                        has_dh = any([hit_id.startswith("PKS_DH") for hit_id in hit_ids])
                        has_er = any([hit_id.startswith("PKS_ER") for hit_id in hit_ids])

                        if has_kr and has_dh and has_er:
                            sequence.append(
                                "polyketide|D"
                            )
                        elif has_kr and has_dh:
                            sequence.append(
                                "polyketide|C"
                            )
                        elif has_kr:
                            sequence.append(
                                "polyketide|B"
                            )
                        else:
                            sequence.append(
                                "polyketide|A"
                            )

                    
                    elif any([hit_id.startswith("Condensation") for hit_id in hit_ids]):
                        sequence.append(
                            "peptide|pubchem|0"
                        )

                    else:
                        # Something else...
                        pass
                    
                if len(sequence) >= 2:
                    yield accession, sequence

def parse_motif(i: int, motif: str) -> str:
    """Parse a motif into a Cypher query.
    
    :param i: Index of motif.
    :type i: int
    :param motif: Motif to parse.
    :type motif: str
    :return: Cypher query.
    :rtype: str
    """
    if match := re.match(r"polyketide\|([A-D])", motif):
        accessory_domains = match.group(1)
        # decoration_type = int(match.group(2))

        return (
            f"(u{i + 1}:Motif {{"
            f"identifier: 'polyketide', "
            f"accessory_domains: '{accessory_domains}', "
            f"decoration_type: '{0}'"
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
    args = cli()
    
    # Configure logger and set log level.
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    # Read all json file paths in the input directory.
    file_paths = glob.glob(f"{args.input}/*.json")
    logger.debug(f"Found {len(file_paths)} files in the input directory.")

    if args.authentication:
        db = GraphDatabase.driver(f"bolt://localhost:{args.port}", auth=tuple(args.authentication))
    else:
        db = GraphDatabase.driver(f"bolt://localhost:{args.port}")


    with db.session() as session:

        # Loop over files and populate the database.
        for file_path in tqdm(file_paths):
            for accession, (identifier, sequence) in enumerate(parse_json(file_path)):

                # Create ProtoCluster node.
                query = """
                CREATE (p:ProtoCluster {identifier: $identifier, accession: $accession})
                """
                session.run(query, identifier=identifier, accession=accession)
                    
                # Create biosynthetic fingerprint node sequence.
                query_begin = """
                CREATE (s:PrimarySequence {identifier: $identifier, accession: $accession, applied_reactions: $applied_reactions})-[:START]->
                """
                query_end = "-[:NEXT]->".join([parse_motif(i, motif) for i, motif in enumerate(sequence)])
                session.run(
                    query_begin + query_end,
                    identifier=identifier,
                    accession=accession,
                    applied_reactions=[]
                )

                # Add link between biosynthetic fingerprint and proto-cluster.
                query = """
                MATCH (p:ProtoCluster {identifier: $identifier})
                MATCH (s:PrimarySequence {identifier: $identifier})
                CREATE (p)-[:HAS_PRIMARY_SEQUENCE]->(s)
                """
                session.run(query, identifier=identifier)


if __name__ == "__main__":
    main()
