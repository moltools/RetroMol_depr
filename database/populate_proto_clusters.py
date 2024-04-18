#!/usr/bin/env python3
"""This script populates the Neo4j database with protoclusters."""
import argparse
from neo4j import GraphDatabase

def cli() -> argparse.Namespace:
    """Parse command-line arguments.
    
    :return: Command-line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to protoclusters.")
    parser.add_argument("--port", default=7687, type=int, help="Neo4j port.")
    parser.add_argument(
        "--authentication",
        default=None,
        nargs=2,
        help="Neo4j authentication as '<username> <password>'."
    )
    return parser.parse_args()

def main() -> None:
    """Driver function."""
    args = cli()
    print(args)

if __name__ == "__main__":
    main()
