# -*- coding: utf-8 -*-

"""Script for investigation speed-up pairwise alignment in database."""

import argparse
import logging
import typing as ty
import time
import json
import re

from neo4j import GraphDatabase
from tqdm import tqdm
from versalign.motif import Motif 
from versalign.sequence import Sequence
from versalign.pairwise import PairwiseAlignment, align_pairwise

from retromol.alignment import PolyketideMotif, PeptideMotif

GLOBAL = PairwiseAlignment.NEEDLEMAN_WUNSCH
LOCAL = PairwiseAlignment.SMITH_WATERMAN


def parse_polyketide_motif(motif: str) -> PolyketideMotif:
    """Parse a polyketide motif.

    :param motif: Motif to parse.
    :type motif: str
    :return: Polyketide motif.
    :rtype: PolyketideMotif
    """
    if match := re.match(r"polyketide\|([A-D])(\d{1,2})", motif):
        accessory_domains = match.group(1)
        decoration_type = int(match.group(2))

        return PolyketideMotif(accessory_domains, decoration_type)

    raise ValueError(f"Invalid polyketide motif: {motif}")


def parse_peptide_motif(motif: str) -> PeptideMotif:
    """Parse a peptide motif.

    :param motif: Motif to parse.
    :type motif: str
    :return: Peptide motif.
    :rtype: PeptideMotif
    """
    if match := re.match(r"peptide\|(\w+)\|(.+)", motif):
        source = match.group(1)
        cid = match.group(2)

        return PeptideMotif(source, cid)

    raise ValueError(f"Invalid peptide motif: {motif}")


def cli() -> argparse.Namespace:
    """Parse command-line arguments.

    :return: Command-line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--port", 
        required=False,
        default=7687, 
        type=int, 
        help="Neo4j port."
    )
    parser.add_argument(
        "--authentication",
        required=True,
        default=None,
        nargs=2,
        help="Neo4j authentication as '<username> <password>'.",
    )
    return parser.parse_args()


# def record_to_sequence(session, record) -> Sequence:
#     identifier = record["b"]["identifier"]

#     result = session.run(
#         (
#             # """
#             # MATCH (:PrimarySequence {identifier: $identifier})-[:START]->(startUnit)
#             # MATCH s = (startUnit)-[:NEXT*]->(endUnit)
#             # WITH nodes(s) AS nodes, length(s) AS path_length
#             # ORDER BY path_length DESC
#             # LIMIT 1
#             # RETURN nodes
#             # """
#             """
#             MATCH (:PrimarySequence {identifier: $identifier})-[:START]->(startMotif)
#             MATCH s = (startMotif)-[:NEXT*]->(endMotif)
#             WHERE NOT (endMotif)-[:NEXT]->()
#             WITH nodes(s) AS motifs
#             RETURN motifs
#             """
#         ), 
#         identifier=identifier
#     )

#     for record in result: # an item could have multiple motif codes
#         motifs = record["motifs"]

    

#     return Sequence("target", [])


def main() -> None:
    """Drive function."""
    args = cli()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")
    logger = logging.getLogger(__name__)
    
    db = GraphDatabase.driver(f"bolt://localhost:{args.port}", auth=tuple(args.authentication))

    start_time = time.time()

    with db.session() as session:

        # result = session.run("""
        #     MATCH (b:PrimarySequence)
        #     WHERE (b)<-[:HAS_PRIMARY_SEQUENCE]-(:Compound) 
        #     OR (b)<-[:HAS_PRIMARY_SEQUENCE]-(:ProtoCluster)
        #     RETURN b""",
        #     fetch_size=1
        # )

        # result = session.run("""
        #     MATCH (b:PrimarySequence)
        #     WHERE (b)<-[:HAS_PRIMARY_SEQUENCE]-(:Compound) 
        #     RETURN b""",
        #     fetch_size=1
        # )

        result = session.run("""
            MATCH (b:TestPrimarySequence)
            RETURN b""",
            fetch_size=1
        )

        query = Sequence("query", [
            PolyketideMotif("B", 7),
            PolyketideMotif("B", 2),
            PolyketideMotif("A", 2),
            PolyketideMotif("D", 7),
            PolyketideMotif("B", 2),
            PolyketideMotif("B", 2),
        ])

        top_scores = []
        top_seqs = []

        for record in tqdm(result):
            # sequence: Sequence = record_to_sequence(session, record)

            name = record["b"]["identifier"]
            
            sequences_str = record["b"]["sequences"]

            sequences_data = json.loads(sequences_str)
            motifs = []
            for datum in sequences_data:
                motif_code = datum["motif_code"]
                for motif in motif_code:
                    if motif.startswith("polyketide"):
                        motifs.append(parse_polyketide_motif(motif))
                    elif motif.startswith("peptide"):
                        motifs.append(parse_peptide_motif(motif))
                    else:
                        raise ValueError(f"Invalid motif: {motif}")
            sequence = Sequence(name, motifs)

            def score_func(a: Motif, b: Motif) -> int:
                if a == b:
                    return 1
                return -1

            options = {
                "gap_penalty": 2,
                # "end_gap_penalty": 1,
            }
            _, _, score = align_pairwise(query, sequence, score_func, LOCAL, options)

            if name == "NPA009987":
                # print(sequences_str)
                print(query)
                print(sequence, score)

            if len(top_scores) < 10:
                top_scores.append(score)
                top_seqs.append(sequence)
            
            elif score > min(top_scores):
                index = top_scores.index(min(top_scores))
                top_scores[index] = score
                top_seqs[index] = sequence
            
            else:
                continue

        # sort
        top_items = sorted(zip(top_seqs, top_scores), key=lambda x: x[1], reverse=True)

        for top_seq, top_score in top_items:
            logger.info(f"{top_seq.identifier} - {top_score} - {top_seq}")

    logger.info(f"Elapsed time: {time.time() - start_time:.2f} s")

    exit(0)


if __name__ == "__main__":
    main()
