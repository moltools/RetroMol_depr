import argparse 
import json 
import os 

from neo4j import GraphDatabase
from tqdm import tqdm

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, required=True, help="Path to dir containing MIBiG JSON files.")
    parser.add_argument("--port", default=7687, type=int, help="Neo4j port.")
    parser.add_argument("--authentication", default=None, nargs=2, help="Neo4j authentication as '<username> <password>'.")
    return parser.parse_args()

def main() -> None:
    args = cli()

    if args.authentication:
        db = GraphDatabase.driver(f"bolt://localhost:{args.port}", auth=tuple(args.authentication))
    else:
        db = GraphDatabase.driver(f"bolt://localhost:{args.port}")

    # Get full paths of all MIBiG JSON files.
    paths = [os.path.join(args.input, f) for f in os.listdir(args.input) if f.endswith(".json")]
    
    with db.session() as session:
        for path in tqdm(paths):

            data = json.load(open(path, "r"))
            data = data["cluster"]

            identifier = data["mibig_accession"]
            organism = data["ncbi_tax_id"]

            # TODO:
            # Get biosyn_class for pathway nodes 
            # Get produced compounds for compound nodes (retrieve pubchem_cid from file)
            # Get bioactivity annotations for bioactivity nodes

            # Get polyketide/peptide motif sequence -- peptide needs to have AA sequence

            print(data.keys())
            
            if "polyketide" in data:
                print(data["polyketide"])

            elif "nrp" in data:
                print(data["nrp"])

            # exit()

    db.close()

if __name__ == "__main__":
    main()
