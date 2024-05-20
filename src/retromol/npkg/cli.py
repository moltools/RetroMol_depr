# -*- coding: utf-8 -*-

"""Command line interface for :mod:`retromol.npkg`."""

import argparse
import logging
import os

from retromol.npkg.connection import Neo4jConnection
from retromol.npkg.purge import purge_database
from retromol.npkg.setup import create_database


def add_subparsers(parser: argparse._SubParsersAction) -> None:
    """Parse command line arguments.

    :param parser: The parser to add subparsers to.
    :type parser: argparse._SubParsersAction
    """
    subparser_create = parser.add_parser("npkg_create")
    subparser_purge = parser.add_parser("npkg_purge")

    for parser in [subparser_create, subparser_purge]:
        parser.add_argument("--uri", type=str, default="bolt://localhost:7687", help="The URI for the Neo4j database.")  # noqa: E501
        parser.add_argument("--usr", type=str, default="neo4j", help="The user for the Neo4j database.")  # noqa: E501
        parser.add_argument("--pwd", type=str, default="password", help="The password for the Neo4j database.")  # noqa: E501

    subparser_create.add_argument("--prc", type=int, default=1, help="Number of processes to use.")  # noqa: E501
    subparser_create.add_argument("--dir", type=str, required=True, help="Path to dir containing data for database construction.")  # noqa: E501
    subparser_create.add_argument("--set-constraints", action="store_true", help="Set constraints on the database.")  # noqa: E501
    subparser_create.add_argument("--only-retrosynthesis", action="store_true", help="Only reparse retrosynthesis data.")  # noqa: E501

    repository_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
    fixtures_path = os.path.join(repository_path, "tests/fixtures")
    path_to_rxn = os.path.join(fixtures_path, "reactions.json")
    path_to_mon = os.path.join(fixtures_path, "monomers.json")

    subparser_create.add_argument("--rxn", type=str, required=False, default=path_to_rxn, help="Path to JSON file containing reaction.")  # noqa: E501
    subparser_create.add_argument("--mon", type=str, required=False, default=path_to_mon, help="Path to JSON file containing monomer patterns.")  # noqa: E501


def main(args: argparse.Namespace) -> None:
    """Driver function.
    
    :param args: The parsed command line arguments.
    :type args: argparse.Namespace
    """
    logger = logging.getLogger(__name__)

    conn = Neo4jConnection(args.uri, args.usr, args.pwd)
    
    if args.mode == "npkg_create":
        create_database(
            conn, 
            path_to_data=args.dir, 
            path_to_rxn=args.rxn,
            path_to_mon=args.mon,
            set_constraints=args.set_constraints,
            only_retrosynthesis=args.only_retrosynthesis,
            num_workers=args.prc
        )

    elif args.mode == "npkg_purge":
        purge_database(conn)

    exit(0)
