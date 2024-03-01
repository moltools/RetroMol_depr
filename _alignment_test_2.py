import typing as ty

from neo4j import GraphDatabase, Session
from tqdm import tqdm 

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

    if record is not None:
        if nodes := record["nodes"]:
            for node in nodes:
                props = dict(node.items())
                seq.append(props)

    return seq 

with db.session() as session:
    print("Connected to Neo4j")

    gap_cost = 2
    end_gap_cost = 1

    # id_a = "NPA009987"
    id_a = "NPA005497"
    seq_a = ModuleSequence(id_a, parse_primary_sequence(retrieve_primary_sequence(session, id_a)))

    # Return 
    query = """
    MATCH (a:PrimarySequence {identifier: $identifier})
    MATCH (b:PrimarySequence)
    WHERE a <> b
    RETURN b
    """
    result = session.run(query, identifier=id_a, fetch_size=1)

    empties = 0

    top_10 = []

    for record in tqdm(result):
        id_b = record["b"]["identifier"]
        seq_b = retrieve_primary_sequence(session, id_b)
        if len(seq_b) != 0:
            seq_b = ModuleSequence(id_b, parse_primary_sequence(seq_b))
            score = seq_a.optimal_alignment(seq_b, gap_cost, end_gap_cost).score

            # Keep 10 best scores.
            if len(top_10) < 10:
                top_10.append((id_b, score))
                top_10.sort(key=lambda x: x[1], reverse=True)
            else:
                if score > top_10[-1][1]:
                    top_10.pop()
                    top_10.append((id_b, score))
                    top_10.sort(key=lambda x: x[1], reverse=True)

        else:
            empties += 1
            continue

    print(f"Empties: {empties}")

    print(top_10)

    # Get associated bioactivities, if any, for the top 10 sequences.
    print("Top 10 sequences and their bioactivities:")
    for id_b, score in top_10:
        query = """
        MATCH (a:PrimarySequence {identifier: $identifier})-[:HAS_BIOACTIVITY]->(b:Bioactivity)
        RETURN b
        """
        result = session.run(query, identifier=id_b)
        for record in result:
            print(record["b"])

db.close()