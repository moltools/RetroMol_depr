# -*- coding: utf-8 -*-

"""Parsing data for the NPKG database."""

import json
import logging
import os
from collections import defaultdict

from rdkit import Chem, RDLogger
from tqdm import tqdm

from retromol.npkg.connection import Neo4jConnection
from retromol.npkg.helpers import safe_get
from retromol.npkg.nodes import (
    BioactivityLabel, 
    BiosyntheticGeneCluster,
    Compound, 
    Motif,
    MotifCode,
    Organism, 
    Pathway
)
from retromol.retrosynthesis.chem import Molecule
from retromol.retrosynthesis.parsing import (
    parse_molecular_patterns,
    parse_reaction_rules,
    parse_mol
)


def parse_compounds(
    conn: Neo4jConnection,
    path_to_rxn: str,
    path_to_mon: str,
) -> None:
    """Parse compounds from the Neo4j database with RetroMol.
    
    :param conn: The Neo4j connection.
    :type conn: Neo4jConnection
    """
    logger = logging.getLogger(__name__)

    # Disable RDKit logger.
    RDLogger.DisableLog("rdApp.*")

    # Parse reactions.
    logger.info("Parsing reaction rules...")
    reactions_src = json.load(open(path_to_rxn, "r", encoding="utf-8"))
    reactions = parse_reaction_rules(reactions_src)
    logger.info(f"Parsed {len(reactions)} reaction rules.")

    # Parse monomers.
    logger.info("Parsing molecular patterns...")
    monomers_src = json.load(open(path_to_mon, "r", encoding="utf-8"))
    monomers = parse_molecular_patterns(monomers_src)
    logger.info(f"Parsed {len(monomers)} molecular patterns.")

    # Retrieve compoounds from 'Polyketides' and 'Amino acids and Peptides' pathways.
    results = conn.query(
        (
            "MATCH (c:Compound)-[:BELONGS_TO]->(p:Pathway) "
            "WHERE p.name IN ['Polyketides', 'Amino acids and Peptides'] "
            "RETURN c.identifier AS identifier, c.source AS source, c.inchi AS inchi"
        ),
        batch_size=1,
    )

    # Unpack results and parse with RetroMol.
    for result in tqdm(results, desc="Parsing compounds with RetroMol", leave=False):
        
        compound_identifier = result["identifier"]
        compound_source = result["source"]
        compound_inchi = result["inchi"]

        # Parse compound with RDKit into RetroMol Molecule.
        try:
            mol = Chem.MolFromInchi(compound_inchi)

            if mol is None:
                raise Exception(f"Failed to parse InChI")

            smi = Chem.MolToSmiles(mol)
            rec = Molecule(f"{compound_source}|{compound_identifier}", smi)

        except Exception as e:
            logger.error(f"Failed to parse compound {compound_identifier} with error: {e}")
            continue

        # Apply retrosynthesis to the compound record.
        # TODO


def parse_donphan(conn: Neo4jConnection, path: str) -> None:
    """Parse the DONPHAN data.
    
    :param conn: The Neo4j connection.
    :type conn: Neo4jConnection
    :param path: The path to the DONPHAN data in csv format.
    :type path: str
    """
    # Disable RDKit logger.
    RDLogger.DisableLog("rdApp.*")

    # Collect bioactivity data.
    bioactivity_library = defaultdict(set)

    with open(path, "r", encoding="utf-8") as file_open:
        header = file_open.readline().strip().split(",")

        for line in tqdm(file_open, desc="Parsing DONPHAN data", leave=False):
            items = list(zip(header, line.strip().split(",")))

            smiles = items[1][1]
            inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))
            inchikey_connectivity = inchikey.split("-")[0]

            bioactivities = items[2:]
            bioactivities = set([h for h, b in bioactivities if b != ""])
            bioactivities = [x[:-3] if x.endswith("_np") else x for x in bioactivities]

            bioactivity_library[inchikey_connectivity].update(bioactivities)

    # Loop through the bioactivity library and add bioactivity nodes and relationships.
    for inchikey_connectivity, bioactivities in tqdm(bioactivity_library.items(), desc="Parsing DONPHAN data", leave=False):  # noqa: E501

        # First retrieve compounds with the given connectivity.
        results = conn.query(
            (
                "MATCH (c:Compound {inchikey_connectivity: $inchikey_connectivity}) "
                "RETURN c.identifier AS identifier, c.source AS source"
            ),
            {"inchikey_connectivity": inchikey_connectivity}
        )
        
        # Annotate retrieved compounds with bioactivity data.
        for result in results:
            
            compound_identifier = result["identifier"]
            compound_source = result["source"]

            for bioactivity in bioactivities:
                bioactivity_source = "donphan"
                
                # Create bioactivity node.
                BioactivityLabel.create(conn, name=bioactivity, source=bioactivity_source)

                # Create relationship between compound and bioactivity.
                conn.query(
                    (
                        "MATCH (c:Compound {identifier: $compound_identifier, source: $compound_source}) "  # noqa: E501
                        "MATCH (b:BioactivityLabel {name: $bioactivity, source: $bioactivity_source}) "  # noqa: E501
                        "MERGE (c)-[:HAS_BIOACTIVITY]->(b)"
                    ),
                    {
                        "compound_identifier": compound_identifier,
                        "compound_source": compound_source,
                        "bioactivity": bioactivity,
                        "bioactivity_source": bioactivity_source,
                    },
                )


def parse_mibig(conn: Neo4jConnection, path: str) -> None:
    """Parse the MIBIG data.
    
    :param conn: The Neo4j connection.
    :type conn: Neo4jConnection
    :param path: The path to the MIBIG data in json format.
    :type path: str
    """
    paths = [
        os.path.join(path, file_path) 
        for file_path in os.listdir(path) 
        if file_path.endswith(".json")
    ]

    for path in tqdm(paths, desc="Parsing MIBIG data", leave=False):
        
        # Get cluster data.
        data = json.load(open(path, "r", encoding="utf-8"))
        cluster = data["cluster"]

        # Retrieve BGC information.
        bgc_identifier = cluster["mibig_accession"]
        bgc_source = "mibig"

        # Create BGC node.
        BiosyntheticGeneCluster.create(conn, identifier=bgc_identifier, source=bgc_source)

        # Parse produced compounds from the cluster, link BGC to compounds.
        if compounds := cluster.get("compounds", None):
            for compound in compounds:

                # Retrieve annotated bioactivities.
                bioactivities = []
                if chem_acts := compound.get("chem_acts", None):
                    for chem_act in chem_acts:
                        bioactivities.append(chem_act["activity"])

                # Retrieve compound NPAtlas identifier, if any, and use it to 
                # link to the compound.
                if database_ids := compound.get("database_id", None):
                    for database_id in database_ids:
                        
                        compound_source, compound_identifier = database_id.split(":")
                        if compound_source == "npatlas":

                            # Link compound to BGC.
                            conn.query(
                                (
                                    "MATCH (c:Compound {identifier: $compound_identifier, source: $compound_source}) "  # noqa: E501
                                    "MATCH (b:BiosyntheticGeneCluster {identifier: $bgc_identifier, source: $bgc_source}) "  # noqa: E501
                                    "MERGE (b)-[:PRODUCES]->(c) "
                                    "MERGE (c)-[:PRODUCED_BY]->(b)"
                                ),
                                {
                                    "compound_identifier": compound_identifier,
                                    "compound_source": compound_source,
                                    "bgc_identifier": bgc_identifier,
                                    "bgc_source": bgc_source,
                                },
                            )

                            # Link compound to bioactivities.
                            for bioactivity in bioactivities:
                                bioactivity_source = "mibig"

                                # Create bioactivity node.
                                BioactivityLabel.create(conn, name=bioactivity, source=bioactivity_source)  # noqa: E501

                                # Create relationship between compound and bioactivity.
                                conn.query(
                                    (
                                        "MATCH (c:Compound {identifier: $compound_identifier, source: $compound_source}) "  # noqa: E501
                                        "MATCH (b:BioactivityLabel {name: $bioactivity, source: $bioactivity_source}) "  # noqa: E501
                                        "MERGE (c)-[:HAS_BIOACTIVITY]->(b)"
                                    ),
                                    {
                                        "compound_identifier": compound_identifier,
                                        "compound_source": compound_source,
                                        "bioactivity": bioactivity,
                                        "bioactivity_source": bioactivity_source,
                                    },
                                )

                                # Create relationship between BGC and bioactivity.
                                conn.query(
                                    (
                                        "MATCH (b:BiosyntheticGeneCluster {identifier: $bgc_identifier, source: $bgc_source}) "  # noqa: E501
                                        "MATCH (b:BioactivityLabel {name: $bioactivity, source: $bioactivity_source}) "  # noqa: E501
                                        "MERGE (b)-[:PRODUCES_COMPOUND_WITH]->(b)"
                                    ),
                                    {
                                        "bgc_identifier": bgc_identifier,
                                        "bgc_source": bgc_source,
                                        "bioactivity": bioactivity,
                                        "bioactivity_source": bioactivity_source,
                                    },
                                )


def parse_npatlas(conn: Neo4jConnection, path: str) -> None:
    """Parse the NPAtlas data.

    :param conn: The Neo4j connection.
    :type conn: Neo4jConnection
    :param path: The path to the NPAtlas data in json format.
    :type path: str
    """
    with open(path, "r") as file:
        data = json.load(file)

    for record in tqdm(data, desc="Parsing NPAtlas data", leave=False):

        # Get producing organism information.
        organism = record["origin_organism"]
        organism_ncbi_id = organism["taxon"]["ncbi_id"]
        organism_type = organism["type"]
        organism_genus = organism["genus"]
        organism_species = organism["species"]

        # Get pathway information.
        pathways = safe_get(record, ["npclassifier", "pathway_results"], [])
        pathways_source = "npclassifier"

        # Get compound information.
        compound_source = "npatlas"
        compound_identifier = record["npaid"]
        compound_inchikey = record["inchikey"]
        compound_inchikey_connectivity = compound_inchikey.split("-")[0]
        compound_inchi = record["inchi"]

        # Create compound node.
        Compound.create(
            conn,
            identifier=compound_identifier,
            source=compound_source,
            inchikey=compound_inchikey,
            inchikey_connectivity=compound_inchikey_connectivity,
            inchi=compound_inchi,
        )

        # Create organism node, and relationship between compound and organism.
        Organism.create(
            conn,
            ncbi_id=organism_ncbi_id,
            type=organism_type,
            genus=organism_genus,
            species=organism_species,
        )

        conn.query(
            (
                "MATCH (c:Compound {identifier: $identifier, source: $source}) "
                "MATCH (o:Organism {ncbi_id: $ncbi_id}) "
                "MERGE (c)-[:PRODUCED_BY]->(o) "
                "MERGE (o)-[:PRODUCES]->(c)"
            ),
            {
                "identifier": compound_identifier,
                "source": compound_source,
                "ncbi_id": organism_ncbi_id,
            },
        )

        # Create pahtway node, and relationship between compound and pathway.
        for pathway in pathways:
            Pathway.create(
                conn,
                name=pathway,
                source=pathways_source,
            )

            conn.query(
                (
                    "MATCH (c:Compound {identifier: $compound_identifier, source: $compound_source}) "  # noqa: E501
                    "MATCH (p:Pathway {name: $pathway_name, source: $pathway_source}) "  # noqa: E501
                    "MERGE (c)-[:BELONGS_TO]->(p) "
                    "MERGE (p)-[:CONTAINS]->(c)"
                ),
                {
                    "compound_identifier": compound_identifier,
                    "compound_source": compound_source,
                    "pathway_name": pathway,
                    "pathway_source": pathways_source,
                },
            )
