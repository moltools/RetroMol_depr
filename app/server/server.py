# -*- coding: utf-8 -*-

"""This module contains the web API for the application."""

from flask import Flask, Response
import neo4j

from retromol.version import get_version

from routes.common import NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD, fail, warning, success
from routes.parsing import blueprint_parse_submission


app = Flask(__name__)
app.register_blueprint(blueprint_parse_submission)


@app.route("/")
def index() -> Response:
    """Return the index page.

    :return: The index page.
    :rtype: Response
    """
    return app.send_static_file("index.html")


@app.route("/api/fetch_server_status", methods=["GET"])
def fetch_server_status() -> Response:
    """API endpoint for checking the server status.

    :return: The response to return.
    :rtype: Response
    """
    return success("Server is up and running!")


@app.route("/api/fetch_database_status", methods=["GET"])
def fetch_database_status() -> Response:
    """API endpoint for checking the database status.

    :return: The response to return.
    :rtype: Response
    """
    try:
        if NEO4J_USER and NEO4J_PASSWORD:
            driver = neo4j.GraphDatabase.driver(
                NEO4J_URI, 
                auth=(NEO4J_USER, NEO4J_PASSWORD)
            )
        else:
            driver = neo4j.GraphDatabase.driver(NEO4J_URI)

        with driver.session() as session:
            result = session.run("MATCH (n) RETURN count(n) AS count")
            count = result.single()["count"]

        driver.close()

        return success(f"Database is up and running with {count} nodes!")

    except Exception as e:
        return fail(f"No database connection: {e}")
    

@app.route("/api/fetch_version", methods=["GET"])
def fetch_version() -> Response:
    """API endpoint for fetching the version of the application.

    :return: The response to return.
    :rtype: Response
    """
    return success("Version fetched successfully!", {"version": get_version()})


def main() -> None:
    """Run the app locally for development."""
    # Main is only run when the file is run directly during local development.
    # When running in Docker, the app is run directly from the Dockerfile.
    app.run(host="localhost", port=4000, debug=True)


if __name__ == "__main__":
    main()
