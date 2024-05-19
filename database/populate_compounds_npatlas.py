#!/usr/bin/env python
"""This script populates the Neo4j database with compounds from the NPAtlas 
database.
"""
import argparse
import json
import typing as ty
from dataclasses import dataclass

from neo4j import GraphDatabase, Session
from tqdm import tqdm


def cli() -> argparse.Namespace:
    """Parse command-line arguments.

    :return: Command-line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to NPAtlas JSON database file.")
    parser.add_argument("--port", default=7687, type=int, help="Neo4j port.")
    parser.add_argument(
        "--authentication",
        default=None,
        nargs=2,
        help="Neo4j authentication as '<username> <password>'.",
    )
    return parser.parse_args()


@dataclass
class CompoundRecord:
    """Dataclass to store compound information."""

    identifier: str
    inchi: str
    inchikey: str
    biosynthetic_pathway: ty.List[str]
    ncbi_ids: ty.List[str]


def add_compound(session: Session, record: CompoundRecord) -> None:
    """
    Add compound to the Neo4j database.

    :param session: Neo4j session.
    :type session: Session
    :param record: Compound record.
    :type record: CompoundRecord
    :return: None
    :rtype: None
    """
    connectivity = record.inchikey.split("-")[0]

    query = """
    CREATE (c:Compound {identifier: $identifier, connectivity: $connectivity, inchi: $inchi, inchikey: $inchikey})
    """
    session.run(
        query,
        identifier=record.identifier,
        connectivity=connectivity,
        inchi=record.inchi,
        inchikey=record.inchikey,
    )

    for pathway in record.biosynthetic_pathway:

        # Create pathway node.
        query = """
        MERGE (b:Pathway {name: $pathway})
        """
        session.run(query, pathway=pathway)

        # Create relationship between compound and pathway.
        query = """
        MATCH (c:Compound {identifier: $identifier})
        MATCH (b:Pathway {name: $pathway})
        MERGE (c)-[:BELONGS_TO]->(b)
        """
        session.run(query, identifier=record.identifier, pathway=pathway)

    for ncbi_id in record.ncbi_ids:

        # Create organism node.
        query = """
        MERGE (o:Organism {ncbi_id: $ncbi_id})
        """
        session.run(query, ncbi_id=ncbi_id)

        # Create relationship between compound and organism.
        query = """
        MATCH (c:Compound {identifier: $identifier})
        MATCH (o:Organism {ncbi_id: $ncbi_id})
        MERGE (c)-[:PRODUCED_BY]->(o)
        """
        session.run(query, identifier=record.identifier, ncbi_id=ncbi_id)


def main() -> None:
    """Driver function."""
    args = cli()
    data = json.load(open(args.input, "r", encoding="utf-8"))

    if args.authentication:
        db = GraphDatabase.driver(f"bolt://localhost:{args.port}", auth=tuple(args.authentication))
    else:
        db = GraphDatabase.driver(f"bolt://localhost:{args.port}")

    with db.session() as session:
        for compound in tqdm(data):

            try:
                ncbi_id = compound["origin_organism"]["taxon"]["ncbi_id"]
                if ncbi_id is not None:
                    ncbi_ids = [ncbi_id]
                else:
                    ncbi_ids = []

            except Exception:
                ncbi_ids = []

            try:
                pathways = compound["npclassifier"]["pathway_results"]
            except Exception:
                pathways = []

            record = CompoundRecord(
                identifier=compound["npaid"],
                inchi=compound["inchi"],
                inchikey=compound["inchikey"],
                biosynthetic_pathway=pathways,
                ncbi_ids=ncbi_ids,
            )

            add_compound(session, record)

    db.close()


if __name__ == "__main__":
    main()
