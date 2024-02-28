import argparse 
import json 
import typing as ty
from collections import defaultdict
from dataclasses import dataclass
from neo4j import Session, GraphDatabase

from rdkit import Chem, RDLogger
from tqdm import tqdm 

RDLogger.DisableLog("rdApp.*")

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to DONPHAN csv file.")
    parser.add_argument("--port", default=7687, type=int, help="Neo4j port.")
    parser.add_argument("--authentication", default=None, nargs=2, help="Neo4j authentication as '<username> <password>'.")
    return parser.parse_args()

@dataclass
class CompoundBioactivityRecord:
    donphan_id: str
    inchikey: str
    bioactivities: ty.List[str]

def main() -> None:
    args = cli()

    # Parse bioactivity data.
    bioactivity_library = defaultdict(set)
    with open(args.input, "r") as file_open:
        header = file_open.readline().strip().split(",")    

        for line in file_open:
            items = list(zip(header, line.strip().split(",")))
            
            smiles = items[1][1]
            inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))

            bioactivities = items[2:]
            bioactivities = set([h.split("_")[0] for h, b in bioactivities if b != ""])
            
            bioactivity_library[inchikey].update(bioactivities)

    # Connect to Neo4j and see if Compound nodes exist in the bioactivity library. 
    # If so, add bioactivity node and relationship to Compound node.
    if args.authentication:
        db = GraphDatabase.driver(f"bolt://localhost:{args.port}", auth=tuple(args.authentication))
    else:
        db = GraphDatabase.driver(f"bolt://localhost:{args.port}")

    with db.session() as session:
        for conn_hash, bioactivities in tqdm(bioactivity_library.items()):
            query = """
            MATCH (c:Compound {inchikey: $inchikey})
            RETURN c.npatlas_id AS npatlas_id
            """
            result = session.run(query, inchikey=conn_hash)
            retrieved = result.single()
            
            if retrieved is not None:
                npatlas_id = retrieved["npatlas_id"]

                if npatlas_id:
                    for bioactivity in bioactivities:
                        query = """
                        MERGE (b:Bioactivity {name: $bioactivity})
                        """
                        session.run(query, bioactivity=bioactivity)

                        query = """
                        MATCH (c:Compound {npatlas_id: $npatlas_id})
                        MATCH (b:Bioactivity {name: $bioactivity})
                        MERGE (c)-[:HAS_BIOACTIVITY]->(b)
                        """
                        session.run(query, npatlas_id=npatlas_id, bioactivity=bioactivity)

if __name__ == "__main__":
    main()