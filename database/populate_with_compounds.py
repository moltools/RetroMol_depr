import argparse 
import json 
import typing as ty
from dataclasses import dataclass
from neo4j import Session, GraphDatabase

from tqdm import tqdm 

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to NPAtlas JSON database file.")
    parser.add_argument("--port", default=7687, type=int, help="Neo4j port.")
    parser.add_argument("--authentication", default=None, nargs=2, help="Neo4j authentication as '<username> <password>'.")
    return parser.parse_args()

@dataclass
class CompoundRecord:
    npatlas_id: str
    inchi: str
    inchikey: str
    biosynthetic_pathway: ty.List[str]
    ncbi_ids: ty.List[str]

def add_compound(session: Session, record: CompoundRecord) -> None:
    query = """
    CREATE (c:Compound {npatlas_id: $npatlas_id, inchi: $inchi, inchikey: $inchikey})
    """
    session.run(query, npatlas_id=record.npatlas_id, inchi=record.inchi, inchikey=record.inchikey)

    for pathway in record.biosynthetic_pathway:

        # Create pathway node.
        query = """
        MERGE (b:Pathway {name: $pathway})
        """
        session.run(query, pathway=pathway)

        # Create relationship between compound and pathway.
        query = """
        MATCH (c:Compound {npatlas_id: $npatlas_id})
        MATCH (b:Pathway {name: $pathway})
        MERGE (c)-[:BELONGS_TO]->(b)
        """
        session.run(query, npatlas_id=record.npatlas_id, pathway=pathway)

    for ncbi_id in record.ncbi_ids:

        # Create organism node.
        query = """
        MERGE (o:Organism {ncbi_id: $ncbi_id})
        """
        session.run(query, ncbi_id=ncbi_id)

        # Create relationship between compound and organism.
        query = """
        MATCH (c:Compound {npatlas_id: $npatlas_id})
        MATCH (o:Organism {ncbi_id: $ncbi_id})
        MERGE (c)-[:PRODUCED_BY]->(o)
        """
        session.run(query, npatlas_id=record.npatlas_id, ncbi_id=ncbi_id)

def main() -> None:
    args = cli()
    data = json.load(open(args.input, "r"))

    if args.authentication:
        db = GraphDatabase.driver(f"bolt://localhost:{args.port}", auth=tuple(args.authentication))
    else:
        db = GraphDatabase.driver(f"bolt://localhost:{args.port}")
    
    with db.session() as session:
        for compound in tqdm(data):
            
            try:
                ncbi_id = compound["origin_organism"]["taxon"]["ncbi_id"]
                if ncbi_id is not None: ncbi_ids = [ncbi_id]
                else: ncbi_ids = []
            except Exception:
                ncbi_ids = []
            
            try:
                pathways = compound["npclassifier"]["pathway_results"]
            except Exception:
                pathways = []

            record = CompoundRecord(
                npatlas_id=compound["npaid"],
                inchi=compound["inchi"],
                inchikey=compound["inchikey"],
                biosynthetic_pathway=pathways,
                ncbi_ids=ncbi_ids
            )

            add_compound(session, record)

if __name__ == "__main__":
    main()