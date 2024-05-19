#!/usr/bin/env python3
"""This script annotates MIBiG clusters in the Neo4j database."""
import argparse
import json
import os

from neo4j import GraphDatabase
from tqdm import tqdm


def cli() -> argparse.Namespace:
    """Parse command-line arguments.

    :return: Command-line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", type=str, required=True, help="Path to dir containing MIBiG JSON files."
    )
    parser.add_argument("--port", default=7687, type=int, help="Neo4j port.")
    parser.add_argument(
        "--authentication",
        default=None,
        nargs=2,
        help="Neo4j authentication as '<username> <password>'.",
    )
    return parser.parse_args()


def main() -> None:
    """Driver function."""
    args = cli()

    if args.authentication:
        db = GraphDatabase.driver(f"bolt://localhost:{args.port}", auth=tuple(args.authentication))
    else:
        db = GraphDatabase.driver(f"bolt://localhost:{args.port}")

    # Get full paths of all MIBiG JSON files.
    paths = [os.path.join(args.input, f) for f in os.listdir(args.input) if f.endswith(".json")]

    with db.session() as session:
        for path in tqdm(paths):

            data = json.load(open(path, "r", encoding="utf-8"))
            cluster = data["cluster"]

            # Add BGC node.
            mibig_id = cluster["mibig_accession"]

            query = """
            MERGE (m:BGC {identifier: $mibig_id})
            """
            session.run(query, mibig_id=mibig_id)

            organisms = [cluster["ncbi_tax_id"]]
            for ncbi_id in organisms:
                query = """
                MATCH (m:BGC {identifier: $mibig_id})
                MERGE (o:Organism {ncbi_id: $ncbi_id})
                MERGE (m)-[:FOUND_IN]->(o)
                """
                session.run(query, mibig_id=mibig_id, ncbi_id=ncbi_id)

            pathways = []
            mibig_pathways = cluster["biosyn_class"]
            if "Polyketide" in mibig_pathways:
                pathways.append("Polyketides")

            if "NRP" in mibig_pathways:
                pathways.append("Amino acids and Peptides")

            for pathway in pathways:
                query = """
                MATCH (m:BGC {identifier: $mibig_id})
                MERGE (p:Pathway {name: $pathway})
                MERGE (m)-[:BELONGS_TO]->(p)
                """
                session.run(query, mibig_id=mibig_id, pathway=pathway)

            # Add compounds with annotated bioactivites. These will be parsed later on into sequences as well.
            if compounds := cluster.get("compounds", None):
                for compound in compounds:

                    bioactivities = []
                    if chem_acts := compound.get("chem_acts", None):
                        for chem_act in chem_acts:
                            bioactivities.append(chem_act["activity"])

                    if database_ids := compound.get("database_id", None):
                        for database_id in database_ids:
                            database, npatlas_id = database_id.split(":")
                            if database == "npatlas":

                                # Link to BGC.
                                query = """
                                MERGE (c:Compound {npatlas_id: $npatlas_id})
                                MERGE (m:BGC {identifier: $mibig_id})
                                MERGE (c)-[:PRODUCED_BY]->(m)
                                """
                                session.run(query, npatlas_id=npatlas_id, mibig_id=mibig_id)

                                # Link to bioactivities.
                                for bioactivity in bioactivities:
                                    query = """
                                    MATCH (c:Compound {npatlas_id: $npatlas_id})
                                    MERGE (b:Bioactivity {name: $bioactivity})
                                    MERGE (c)-[:HAS_BIOACTIVITY]->(b)
                                    """
                                    session.run(
                                        query, npatlas_id=npatlas_id, bioactivity=bioactivity
                                    )

    db.close()


if __name__ == "__main__":
    main()
