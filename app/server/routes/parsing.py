# -*- coding: utf-8 -*-

"""This module contains functions for parsing input data."""

import json
import os

import neo4j
from flask import Blueprint, Response, request

from retromol.npkg.antismash import get_antismash_data, parse_antismash_json
from retromol.retrosynthesis.parsing import (
    Molecule,
    parse_mol, 
    parse_molecular_patterns, 
    parse_reaction_rules
)

from .common import fail, warning, success


try:
    absolute_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    monomers = json.load(open(os.path.join(absolute_path, "data/monomers.json"), "r", encoding="utf-8"))
    reactions = json.load(open(os.path.join(absolute_path, "data/reactions.json"), "r", encoding="utf-8"))
    REACTIONS = parse_reaction_rules(reactions)
    MONOMERS = parse_molecular_patterns(monomers)
except Exception:
    REACTIONS = []
    MONOMERS = []


blueprint_parse_submission = Blueprint("parse_submission", __name__)
@blueprint_parse_submission.route("/api/parse_submission", methods=["POST"])
def parse_submission() -> Response:
    """Parse submission data."""
    
    # Unpack data.
    data = request.json

    try:
        input_type = data["inputType"]
        input_value = data["inputValue"]

        assert input_type in ["smiles", "jobId", "json"], "Invalid input type."
        assert len(input_value) > 0, "Input value is empty."

    except Exception as e:
        error_type = e.__class__.__name__
        return fail(f"Invalid input data ({error_type}): {e}")
    
    # Parse data.
    try:

        ########################################################################
        #
        # Parse SMILES
        #
        ########################################################################
        if input_type == "smiles":
            molecule = Molecule("input", input_value)
            result = parse_mol(molecule, REACTIONS, MONOMERS)
            sequences = [sequence["motif_code"] for sequence in result.sequences]
            payload = {"sequences": sequences}
            
            if result.success is True:
                return success("Molecule parsed successfully!", payload)
            else:
                return fail("Failed to parse molecule!")
        
        ########################################################################
        #
        # Parse job ID
        #
        ########################################################################
        elif input_type == "jobId":
            data = get_antismash_data(input_value)
            sequences = parse_antismash_json(data)
            payload = {"sequences": sequences}

            return success("Job parsed successfully!", payload)
        
        ########################################################################
        #
        # Parse JSON
        #
        ########################################################################
        elif input_type == "json":
            data = json.loads(input_value)
            sequences = parse_antismash_json(data)
            payload = {"sequences": sequences}

            return success("JSON parsed successfully!")
    
    except Exception as e:
        error_type = e.__class__.__name__
        return fail(f"Failed to parse input ({error_type}): {e}")

    return warning("Endpoint not implemented yet!")
