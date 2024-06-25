import argparse
import json
import numpy as np
import os
from tqdm import tqdm

from retromol.npkg.connection import Neo4jConnection
from versalign.motif import Motif
from versalign.pairwise import PairwiseAlignment, align_pairwise
from versalign.sequence import Sequence


class A(Motif):
    def __eq__(self, other): return isinstance(other, A)
    def __str__(self): return "A"

class B(Motif):
    def __eq__(self, other): return isinstance(other, B)
    def __str__(self): return "B"

class C(Motif):
    def __eq__(self, other): return isinstance(other, C)
    def __str__(self): return "C"

class D(Motif):
    def __eq__(self, other): return isinstance(other, D)
    def __str__(self): return "D"

class Pep(Motif):
    def __eq__(self, other): return isinstance(other, Pep)
    def __str__(self): return "Pep"


def score_func(a, b):
    if a == b:
        return 3
    return -1


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=str, default=".", help="The output directory.")
    parser.add_argument("--uri", type=str, default="bolt://localhost:7687", help="The URI for the Neo4j database.")
    parser.add_argument("--usr", type=str, default="neo4j", help="The user for the Neo4j database.")
    parser.add_argument("--pwd", type=str, default="password", help="The password for the Neo4j database.")
    return parser.parse_args()


def parse_sequence(identifier, seq_src):
    data = json.loads(seq_src)
    seq = []
    for x in data:
        if x.startswith("peptide"): seq.append(Pep())
        elif x.startswith("polyketide"):
            _, motif_code = x.split("|")
            if motif_code.startswith("A"): seq.append(A())
            elif motif_code.startswith("B"): seq.append(B())
            elif motif_code.startswith("C"): seq.append(C())
            elif motif_code.startswith("D"): seq.append(D())
            else: raise ValueError(f"Unknown motif code: {motif_code}")
        else: raise ValueError(f"Unknown motif code: {x}")
    return Sequence(identifier, seq)


def main():
    args = cli()
    conn = Neo4jConnection(args.uri, args.usr, args.pwd)
    labels, identifiers, seqs = [], [], []
    result = conn.query("MATCH (n:MotifCode) RETURN n")
    for record in result:
        label = record["n"]["compound_source"]  # npatlas / asdb4
        identifier = record["n"]["compound_identifier"]
        seq_src = record["n"]["src"]  # jsonified str sequence
        labels.append(label)
        identifiers.append(identifier)
        seqs.append(parse_sequence(identifier, seq_src))

    # Separate labels, identifiers, seqs based on label.
    npatlas_labels, npatlas_identifiers, npatlas_seqs = [], [], []
    asdb4_labels, asdb4_identifiers, asdb4_seqs = [], [], []

    for label, identifier, seq in zip(labels, identifiers, seqs):
        if label == "npatlas":
            npatlas_labels.append(label)
            npatlas_identifiers.append(identifier)
            npatlas_seqs.append(seq)
        elif label == "asdb4":
            asdb4_labels.append(label)
            asdb4_identifiers.append(identifier)
            asdb4_seqs.append(seq)
        else:
            raise ValueError(f"Unknown label: {label}")
        

    query = Sequence("thermolide", [B(), B(), B(), D(), B(), D(), B(), B(), Pep()]) # pep = 602 (alanine)



    exit()

    # Align query sequence with all sequences in the database.
    for j, query in enumerate(asdb4_seqs):
        scores = []
        for seq in tqdm(npatlas_seqs):
            _, _, score = align_pairwise(
                query, 
                seq, 
                score_func, 
                PairwiseAlignment.SMITH_WATERMAN,
                # options={"endGapPenalty": -1},
            )
            scores.append(score)
        top_n = np.argsort(scores)[::-1][:10]

        print(f"Query: {asdb4_identifiers[j]}")
        for i in top_n: print(f"{npatlas_labels[i]}: {npatlas_identifiers[i]}")


if __name__ == "__main__":
    main()
