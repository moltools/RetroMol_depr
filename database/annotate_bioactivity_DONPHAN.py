import argparse 
from collections import defaultdict
from neo4j import GraphDatabase

from rdkit import Chem, RDLogger
from tqdm import tqdm 

RDLogger.DisableLog("rdApp.*")

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to DONPHAN csv file.")
    parser.add_argument("--port", default=7687, type=int, help="Neo4j port.")
    parser.add_argument("--authentication", default=None, nargs=2, help="Neo4j authentication as '<username> <password>'.")
    return parser.parse_args()

def main() -> None:
    args = cli()

    # Parse bioactivity data.
    bioactivity_library = defaultdict(set)
    with open(args.input, "r") as file_open:
        header = file_open.readline().strip().split(",")    

        for line in tqdm(file_open):
            items = list(zip(header, line.strip().split(",")))
            
            smiles = items[1][1]
            inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))
            connectivity = inchikey.split("-")[0]

            bioactivities = items[2:]
            bioactivities = set([h for h, b in bioactivities if b != ""])
            bioactivities = [x[:-3] if x.endswith("_np") else x for x in bioactivities]
            
            bioactivity_library[connectivity].update(bioactivities)

    # Connect to Neo4j and see if Compound nodes exist in the bioactivity library. 
    # If so, add bioactivity node and relationship to Compound node.
    if args.authentication:
        db = GraphDatabase.driver(f"bolt://localhost:{args.port}", auth=tuple(args.authentication))
    else:
        db = GraphDatabase.driver(f"bolt://localhost:{args.port}")

    with db.session() as session:
        for connectivity, bioactivities in tqdm(bioactivity_library.items()):
            query = """
            MATCH (c:Compound {connectivity: $connectivity})
            RETURN c.identifier AS identifier
            """
            result = session.run(query, connectivity=connectivity)
            
            for retrieved in result:
                identifier = retrieved["identifier"]

                if identifier is not None:
                    for bioactivity in bioactivities:
                        query = """
                        MERGE (b:Bioactivity {name: $bioactivity, source: "DONPHAN"})
                        """
                        session.run(query, bioactivity=bioactivity)

                        query = """
                        MATCH (c:Compound {identifier: $identifier})
                        MATCH (b:Bioactivity {name: $bioactivity})
                        MERGE (c)-[:HAS_BIOACTIVITY]->(b)
                        """
                        session.run(query, identifier=identifier, bioactivity=bioactivity)

    db.close()

if __name__ == "__main__":
    main()