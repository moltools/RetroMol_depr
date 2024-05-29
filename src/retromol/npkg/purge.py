# -*- coding: utf-8 -*-

"""Functions for purging the Neo4j database."""

import logging

from retromol.npkg.connection import Neo4jConnection


def remove_constraints(conn: Neo4jConnection) -> None:
    """Remove constraints from the Neo4j database.

    :param conn: The Neo4j connection.
    :type conn: Neo4jConnection
    :raises TypeError: If the connection is not a Neo4jConnection.
    """
    logger = logging.getLogger(__name__)
    logger.info("Removing constraints...")

    if not isinstance(conn, Neo4jConnection):
        msg = f"Expected a Neo4jConnection but received {type(conn)}"
        logger.error(msg)
        raise TypeError(msg)

    # Remove constraints.
    for constraint in conn.query("SHOW constraints"):
        constraint_name = constraint["name"]
        conn.query(f"DROP CONSTRAINT {constraint_name}")

    logger.info("Constraints removed.")


def purge_database(conn: Neo4jConnection) -> None:
    """Purge the Neo4j database.

    :param conn: The Neo4j connection.
    :type conn: Neo4jConnection
    :raises TypeError: If the connection is not a Neo4jConnection.
    """
    logger = logging.getLogger(__name__)
    logger.info("Purging database...")

    if not isinstance(conn, Neo4jConnection):
        msg = f"Expected a Neo4jConnection but received {type(conn)}"
        logger.error(msg)
        raise TypeError(msg)

    # Remove contraints.
    remove_constraints(conn)

    # Purge the database.
    while True:
        query = (
            "MATCH (n) WITH n LIMIT $limit " "DETACH DELETE n " "RETURN count(n) AS deleted_count"
        )
        result = conn.query(query, {"limit": 1000})
        deleted_count = result[0]["deleted_count"]
        if deleted_count == 0:
            break

    logger.info("Database purged.")


def purge_retrosynthesis_data(conn: Neo4jConnection, only_calculated: bool = True) -> None:
    """Purge the retrosynthesis data from the Neo4j database.

    :param conn: The Neo4j connection.
    :type conn: Neo4jConnection
    :param only_calculated: If True, only calculated retrosynthesis data will be purged.
    :type only_calculated: bool
    :raises TypeError: If the connection is not a Neo4jConnection.
    """
    logger = logging.getLogger(__name__)
    logger.info("Purging retrosynthesis data...")

    if not isinstance(conn, Neo4jConnection):
        msg = f"Expected a Neo4jConnection but received {type(conn)}"
        logger.error(msg)
        raise TypeError(msg)

    # Remove contraints.
    remove_constraints(conn)

    # Construct the query for purging database from MotifCodes.
    if only_calculated:
        query = (
            "MATCH (n:MotifCode {calculated: true}) "
            "WITH n LIMIT $limit "
            "DETACH DELETE n "
            "RETURN count(n) AS deleted_count "
            "UNION ALL "
            "MATCH (m:Motif {calculated: true}) "
            "WITH m LIMIT $limit "
            "DETACH DELETE m "
            "RETURN count(m) AS deleted_count "
        )
    else:
        query = (
            "MATCH (n:MotifCode) "
            "WITH n LIMIT $limit "
            "DETACH DELETE n "
            "RETURN count(n) AS deleted_count "
            "UNION ALL "
            "MATCH (m:MotifCode) "
            "WITH m LIMIT $limit "
            "DETACH DELETE m "
            "RETURN count(m) AS deleted_count"
        )

    # Purge the MotifCode retrosynthesis data.
    while True:
        result = conn.query(query, {"limit": 1000})
        deleted_count = result[0]["deleted_count"]
        if deleted_count == 0:
            break

    logger.info("Retrosynthesis data purged.")
