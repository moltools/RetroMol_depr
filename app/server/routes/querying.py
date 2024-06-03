# -*- coding: utf-8 -*-

"""This module contains functions for quering a query."""

import json
import os
import re
import typing as ty

import neo4j
from flask import Blueprint, Response, request
from versalign.motif import Motif
from versalign.pairwise import PairwiseAlignment, align_pairwise
from versalign.sequence import Sequence
from tqdm import tqdm

from retromol.npkg.antismash import get_antismash_data, parse_antismash_json
from retromol.retrosynthesis.alignment import OtherMotif, PolyketideMotif, sequence_from_motif_string_list
from retromol.retrosynthesis.parsing import (
    Molecule,
    parse_mol, 
    parse_molecular_patterns, 
    parse_reaction_rules
)

from .common import NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER, fail, success, warning


def pairwise_match(
    query: ty.List[ty.Dict[str, ty.Any]], 
    gap_penalty: int, 
    end_gap_penalty: int
) -> ty.Dict[str, ty.Any]:
    """Pairwise match query against database.
    
    :param query: Query to match.
    :type query: ty.List[ty.Dict[str, ty.Any]]
    :param gap_penalty: Gap penalty.
    :type gap_penalty: int
    :param end_gap_penalty: End gap penalty.
    :type end_gap_penalty: int
    :return: Match results.
    :rtype: ty.Dict[str, ty.Any]
    """
    if NEO4J_USER and NEO4J_PASSWORD:
        driver = neo4j.GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))
    else:
        driver = neo4j.GraphDatabase.driver(NEO4J_URI)

    top_scores, top_seqs = [], []
    with driver.session() as session:
        result = session.run("MATCH (b:MotifCode) RETURN b", fetch_size=1)
        
        query = Sequence(
            "query",
            [
                # OtherMotif(cid="594"),
                PolyketideMotif("B", 7),
                PolyketideMotif("B", 2),
                PolyketideMotif("A", 2),
                PolyketideMotif("D", 7),
                PolyketideMotif("B", 2),
                PolyketideMotif("B", 2),
            ],
        )

        for record in tqdm(result):
            name = record["b"]["compound_identifier"]
            motif_code = record["b"]["src"]
            motif_code = json.loads(motif_code)
            sequence = sequence_from_motif_string_list(name, motif_code)

            def score_func(a: Motif, b: Motif) -> int:
                if a == b:
                    return 1
                return -1

            max_num_matches = 10
            alignment_type = PairwiseAlignment.SMITH_WATERMAN  # local
            # alignment_type = PairwiseAlignment.NEEDLEMAN_WUNSCH  # global
            options = {"gap_penalty": 2}
            # options = {"gap_penalty": 2, "end_gap_penalty": 1}

            _, _, score = align_pairwise(
                query, 
                sequence, 
                score_func, 
                alignment_type, 
                options
            )

            if len(top_scores) < max_num_matches:
                top_scores.append(score)
                top_seqs.append(sequence)

            elif score > min(top_scores):
                index = top_scores.index(min(top_scores))
                top_scores[index] = score
                top_seqs[index] = sequence

            else:
                continue

    # Sort
    top_items = sorted(zip(top_seqs, top_scores), key=lambda x: x[1], reverse=True)
    
    for top_seq, top_score in top_items:
        print(f"{top_seq.identifier} - {top_score} - {top_seq}")

    return {"top_items": []}


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
        gap_penalty = data["gapPenalty"]
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
    
    if query_type == "match":
        results = pairwise_match(
            query=query,
            gap_penalty=gap_penalty,
            end_gap_penalty=end_gap_penalty,
        )
        return success("Query ran successfully.", results)
    
    elif query_type == "query":
        return warning("Endpoint not implemented yet!")
    
    return fail("Endpoint not implemented yet!")
