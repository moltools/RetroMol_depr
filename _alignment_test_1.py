import typing as ty

from neo4j import GraphDatabase, Session

from retromol_sequencing.alignment import ModuleSequence, parse_primary_sequence, MultipleSequenceAlignment

db = GraphDatabase.driver("bolt://localhost:7687", auth=("neo4j", "password"))

def retrieve_primary_sequence(session: Session, identifier: str) -> ty.List[ty.Dict[str, ty.Any]]:
    query = """
    MATCH (:PrimarySequence {identifier: $identifier})-[:START]->(startUnit)
    MATCH s = (startUnit)-[:NEXT*]->(endUnit)
    WITH nodes(s) AS nodes, length(s) AS path_length
    ORDER BY path_length DESC
    LIMIT 1
    RETURN nodes
    """
    result = session.run(query, identifier=identifier)
    record = result.single()
    
    seq = []
    for node in record["nodes"]:
        props = dict(node.items())
        seq.append(props)

    return seq 

with db.session() as session:
    print("Connected to Neo4j")

    npatlas_id_a = "NPA009987"
    npatlas_id_b = "NPA005497"

    seq_a = retrieve_primary_sequence(session, npatlas_id_a)
    seq_b = retrieve_primary_sequence(session, npatlas_id_b)

    seq_a = ModuleSequence(npatlas_id_a, parse_primary_sequence(seq_a))
    seq_b = ModuleSequence(npatlas_id_b, parse_primary_sequence(seq_b))

    gap_cost = 2
    end_gap_cost = 1
    alignment = seq_a.optimal_alignment(seq_b, gap_cost, end_gap_cost)

    print(alignment.score)
    print(alignment.seq1)
    print(alignment.seq2)

    print(len(alignment.seq1))
    print(len(alignment.seq2))

    print(seq_a.optimal_alignment(seq_a, gap_cost, end_gap_cost).score)

    msa = MultipleSequenceAlignment([seq_a, seq_b, seq_a], gap_cost, end_gap_cost).get_alignment()
    print(msa)

db.close()