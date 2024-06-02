# -*- coding: utf-8 -*-

"""This module contains functions for quering a query."""

import json
import os
import re
import typing as ty

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


# def motif_code_to_query(motif_code: ty.List[str]) -> ty.List[ty.List[str]]:
#     """Convert motif code to a Cypher query.
    
#     :param motif_code: List of motifs.
#     :type motif_code: list[str]
#     :return: List of queries.
#     :rtype: list[list[str]]
#     """
#     wrapped_motif_code = []

#     for motif in motif_code:
#         if match := re.match(r"polyketide\|([A-D])(\d{1,2})", motif):
#             polyketide_type = match.group(1)
#             polyketide_decor = int(match.group(2))
            
#             parsed_motif = {
#                 "motifType": "polyketide",
#                 "polyketideType": polyketide_type,
#                 "polyketideDecor": polyketide_decor,
#                 "peptideSource": "Any",
#                 "peptideCid": "Any"
#             }

#         elif match := re.match(r"peptide\|(\w+)\|(.+)", motif):
#             peptide_source = match.group(1)
#             peptide_cid = int(match.group(2))

#             parsed_motif = {
#                 "motifType": "peptide",
#                 "polyketideType": "Any",
#                 "polyketideDecor": "Any",
#                 "peptideSource": peptide_source,
#                 "peptideCid": peptide_cid
#             }

#         else:
#             raise ValueError(f"Invalid motif: {motif}")    

#         wrapped_motif_code.append([parsed_motif])  # Wrap individual motifs in a list.
    
#     return wrapped_motif_code


blueprint_query_submission = Blueprint("query_submission", __name__)
@blueprint_query_submission.route("/api/query_submission", methods=["POST"])
def query_submission() -> Response:
    """Query submission data."""
    
    # Unpack data.
    data = request.json

    try:
        query = data["query"]
        query_type = data["queryType"]
        alignment_type = data["alignmentType"]
        end_gap_penalty = data["endGapPenalty"]
        gap_enalty = data["gapPenalty"]
        min_match_length = data["minMatchLength"]
        max_match_length = data["maxMatchLength"]
        max_num_matches = data["maxNumMatches"]
        query_against_molecules = data["queryAgainstMolecules"]
        query_against_protoclusters = data["queryAgainstProtoclusters"]
        query_has_leading_modules = data["queryHasLeadingModules"]
        query_has_trailing_modules = data["queryHasTrailingModules"]
        selected_bioactivity_labels = data["selectedBioactivityLabels"]
        selected_organism_labels = data["selectedOrganismLabels"]

        assert len(query) > 0, "Input query is empty."
        assert query_type in ["match", "query"], "Invalid query type."
        assert alignment_type in ["local", "global"], "Invalid alignment strategy."
        assert min_match_length > 0, "Invalid minimum match length."
        assert max_match_length > 0, "Invalid maximum match length."
        assert max_match_length >= min_match_length, "Invalid match length range."
        assert max_num_matches > 0, "Invalid maximum number of matches."
        assert max_num_matches <= 100, "Maximum number of matches is too high."

    except Exception as e:
        error_type = e.__class__.__name__
        return fail(f"Invalid input data ({error_type}): {e}")
    
    return success("Endpoint not implemented yet!")
