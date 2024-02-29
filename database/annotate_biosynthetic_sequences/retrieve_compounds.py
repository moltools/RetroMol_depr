import argparse 

from neo4j import GraphDatabase
from rdkit import Chem, RDLogger
from tqdm import tqdm

RDLogger.DisableLog("rdApp.*")

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rules", type=str, required=True, help="Path to JSON file containing rules.")
    parser.add_argument("--out", type=str, required=True, help="Path to output file.")
    parser.add_argument("--port", default=7687, type=int, help="Neo4j port.")
    parser.add_argument("--authentication", default=None, nargs=2, help="Neo4j authentication as '<username> <password>'.")
    return parser.parse_args()

def main() -> None:
    args = cli()
    
    if args.authentication:
        db = GraphDatabase.driver(f"bolt://localhost:{args.port}", auth=tuple(args.authentication))
    else:
        db = GraphDatabase.driver(f"bolt://localhost:{args.port}")

    out_file = open(args.out, "w")
    out_file.write("npatlas_id\tsmiles\n")

    with db.session() as session:
        query = """
        MATCH (c:Compound)-[:BELONGS_TO]->(p:Pathway {name: 'Polyketides'})
        RETURN c.npatlas_id AS npatlas_id, c.inchi AS inchi
        """
        result = session.run(query, fetch_size=100)

        for record in tqdm(result):
            npatlas_id = record["npatlas_id"]
            inchi = record["inchi"]

            try:
                smiles = Chem.MolToSmiles(Chem.MolFromInchi(inchi))
            except Exception:
                continue
            
            out_file.write(f"{npatlas_id}\t{smiles}\n")

    out_file.close()

if __name__ == "__main__":
    main()