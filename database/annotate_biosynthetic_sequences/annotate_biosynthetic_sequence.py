import argparse 
import os 
import json
import typing as ty
from neo4j import GraphDatabase

from tqdm import tqdm 
from rdkit import RDLogger

from retromol.parsing import Result
from retromol_sequencing.sequencing import parse_modular_natural_product

RDLogger.DisableLog("rdApp.*")

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=str, required=True, help="Path to dir containing JSON results.")
    parser.add_argument("--port", default=7687, type=int, help="Neo4j port.")
    parser.add_argument("--authentication", default=None, nargs=2, help="Neo4j authentication as '<username> <password>'.")
    return parser.parse_args()

def parse_motif(i: int, motif: ty.Dict[str, ty.Any]) -> str:

    motif_type = motif["identifier"]

    if motif_type == "polyketide":
        props = motif["properties"]
        accessory_domains = props["accessory_domains"]
        decoration_type = props["decoration_type"]

        return f"(u{i + 1}:Motif {{identifier: '{motif_type}', accessory_domains: {accessory_domains}, decoration_type: '{decoration_type}'}})"
    
    elif motif_type == "peptide":
        props = motif["properties"]
        pubchem_cid = props["pubchem_cid"]
        classification = props["classification"]

        return f"(u{i + 1}:Motif {{identifier: '{motif_type}', pubchem_cid: '{pubchem_cid}', classification: '{classification}'}})"
    
    else:
        return f"(u{i + 1}:Motif {{identifier: '{motif_type}'}})"

def main() -> None:
    args = cli()

    # Get all paths to JSON files in results directory.
    paths = [os.path.join(args.results, path) for path in os.listdir(args.results) if path.endswith(".json")]
    print(f"Found {len(paths)} JSON files in {args.results}.")

    if args.authentication:
        db = GraphDatabase.driver(f"bolt://localhost:{args.port}", auth=tuple(args.authentication))
    else:
        db = GraphDatabase.driver(f"bolt://localhost:{args.port}")

    with db.session() as session:
        
        # Iterate over all JSON files.
        for path in tqdm(paths):
            
            result = Result.from_json(path)
        
            if result.success == False:
                continue
            
            if not result.has_identified_monomers():
                continue

            identifier = result.identifier # We assume that the identifier is the NPAtlas ID.
            applied_reactions = result.applied_reactions
            sequences = parse_modular_natural_product(result)
            
            # Add biosynthetic fingerprint to graph.
            for accession, sequence in enumerate(sequences):

                # Create biosynthetic fingerprint node sequence.
                query_begin = """
                CREATE (s:PrimarySequence {identifier: $identifier, accession: $accession, applied_reactions: $applied_reactions})-[:START]->
                """
                query_end = "-[:NEXT]->".join([parse_motif(i, motif) for i, motif in enumerate(sequence)])
                session.run(query_begin + query_end, identifier=identifier, accession=accession, applied_reactions=applied_reactions)

                # Add link between biosynthetic fingerprint and compound.
                query = """
                MATCH (c:Compound {identifier: $identifier})
                MATCH (s:PrimarySequence {identifier: $identifier})
                CREATE (c)-[:HAS_PRIMARY_SEQUENCE]->(s)
                """
                session.run(query, identifier=identifier)

    db.close()

if __name__ == "__main__":
    main()