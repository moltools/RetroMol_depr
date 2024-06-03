# -*- coding: utf-8 -*-

"""Functions for parsing AntiSMASH output."""

import os
import logging
import re
import typing as ty
from collections import defaultdict

from Bio.Seq import Seq
import requests

from retromol.npkg.paras import label_to_pubchem_cid, predict_specificity


class Gene:
    """Gene in AntiSMASH data."""

    def __init__(
        self, 
        name: str,
        start: int, 
        end: int, 
        strand: str, 
    ) -> None:
        """Initialize a gene.
        
        :param name: Name of the gene.
        :type name: str
        :param start: Start of the gene.
        :type start: int
        :param end: End of the gene.
        :type end: int
        :param strand: Strand of the gene.
        :type strand: str
        """
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.modules = list()

    def __str__(self) -> str:
        """Get the string representation of the gene.
        
        :return: String representation of the gene.
        :rtype: str
        """
        return f"Gene(name={self.name}, start={self.start}, end={self.end}, strand={self.strand})"

    def __hash__(self) -> int:
        """Get the hash of the locus.
        
        :return: Hash of the locus.
        :rtype: int
        """
        return hash(self.name)
    
    def __eq__(self, other: ty.Any) -> bool:
        """Check if two genes are equal.

        :param other: Other gene.
        :type other: ty.Any
        :return: True if the genes are equal, False otherwise.
        :rtype: bool
        """
        if not isinstance(other, Gene):
            return False
        return self.name == other.name


class ProtoCluster:
    """ProtoCluster in AntiSMASH data."""

    def __init__(
        self,
        start: int,
        end: int,
        product: str,
        category: str,
    ) -> None:
        """Initialize a proto cluster.
        
        :param start: Start of the proto cluster.
        :type start: int
        :param end: End of the proto cluster.
        :type end: int
        :param product: Product of the proto cluster.
        :type product: str
        :param category: Category of the proto cluster.
        :type category: str
        """
        self.start = start
        self.end = end
        self.product = product
        self.category = category
        self.genes = list()

    def __str__(self) -> str:
        """Get the string representation of the proto cluster.
        
        :return: String representation of the proto cluster.
        :rtype: str
        """
        return f"ProtoCluster(start={self.start}, end={self.end}, product={self.product}, category={self.category})"

    def __hash__(self) -> int:
        """Get the hash of the proto cluster.
        
        :return: Hash of the proto cluster.
        :rtype: int
        """
        return hash((self.start, self.end, self.product, self.category))
    
    def __eq__(self, other: ty.Any) -> bool:
        """Check if two proto clusters are equal.
        
        :param other: Other proto cluster.
        :type other: ty.Any
        :return: True if the proto clusters are equal, False otherwise.
        :rtype: bool
        """
        if not isinstance(other, ProtoCluster):
            return False
        return (
            self.start == other.start and
            self.end == other.end and
            self.product == other.product and
            self.category == other.category
        )
    
    def __len__(self) -> int:
        """Get the number of genes in the proto cluster.
        
        :return: Number of genes in the proto cluster.
        :rtype: int
        """
        return len(self.genes)
    
    def to_motif_code(self, dna: str) -> ty.Dict[str, ty.Any]:
        """Convert the proto cluster to a motif code.
        
        :param dna: DNA sequence of source.
        :type dna: str
        :return: Motif code.
        :rtype: ty.Dict[str, ty.Any]
        """
        logger = logging.getLogger(__name__)

        # Collect sequences for PARAS A-domain specificity prediction.
        to_predict = []

        # Bundle relevant modules into putative domains.
        groups = []
        for gene in self.genes:
            gene_sequence_added = False
            modules = gene.modules

            current_group = []
            for module in modules:
                if module in ["PKS_KS", "AMP-binding"]:
                    if module == "AMP-binding" and not gene_sequence_added:
                        gene_start = gene.start
                        gene_end = gene.end
                        gene_strand = gene.strand
                        gene_dna = dna[gene_start:gene_end]
                        if gene_strand == "-":
                            gene_dna = str(Seq(gene_dna).reverse_complement())
                        gene_protein = str(Seq(gene_dna).translate())
                        gene_fasta = f">{gene.name}\n{gene_protein}"
                        to_predict.append(gene_fasta)
                        gene_sequence_added = True

                    if current_group:
                        groups.append(current_group)

                    current_group = [module]

                elif (
                    module in ["PKS_DH", "PKS_ER", "PKS_KR"] 
                    and current_group 
                    and current_group[0] == "PKS_KS"
                ):
                    current_group.append(module)

                else:
                    continue 

            if current_group:
                groups.append(current_group)

        # Predict specificity for A-domains.
        to_predict = "\n".join(to_predict)
        try:
            predicted_specificities = predict_specificity(to_predict)
            predicted_specificities = [p[1][0][1] for p in predicted_specificities]
        except Exception as e:
            logger.error(f"Could not predict A-domain specificities: {e}")
            predicted_specificities = []

        # Count number of groups with AMP-domains.
        num_amp_binding = len([group for group in groups if "AMP-binding" in group])
        if num_amp_binding != len(predicted_specificities):
            predicted_specificities = ["Any"] * num_amp_binding

        predicted_specificities = iter(predicted_specificities)
        
        # Parse the putative domains.
        query = []
        for group in groups:
            if "PKS_KS" in group and "AMP-binding" in group:
                raise ValueError("Found both PKS_KS and AMP-binding in the same domain.")
            
            elif "PKS_KS" in group:
                polyketide_type = "A"
                polyketide_decor = "Any"

                if "PKS_ER" in group:
                    polyketide_type = "D"
                
                elif "PKS_DH" in group:
                    polyketide_type = "C" 
                
                elif "PKS_KR" in group:
                    polyketide_type = "B"  
                
                query.append([{
                    "motifType": "polyketide",
                    "polyketideType": polyketide_type,
                    "polyketideDecor": polyketide_decor,
                    "peptideSource": "Any",
                    "peptideCid": "Any"
                }])

            elif "AMP-binding" in group:
                specificity = next(predicted_specificities)
                pubchem_cid = label_to_pubchem_cid(specificity)
                if pubchem_cid is None:
                    pubchem_cid = "Any"

                query.append([{
                    "motifType": "peptide",
                    "polyketideType": "Any",
                    "polyketideDecor": "Any",
                    "peptideSource": "Any",
                    "peptideCid": pubchem_cid
                }])

            else: 
                raise ValueError("Could not find PKS_KS or AMP-binding in the domain.")

        return {   
            "title": f"{self.category} ({self.product})",
            "queryType": "parsed",
            "query": query,
            "metadata": {}
        }


class ProtoClusterDict:
    """Proto cluster dictionary for AntiSMASH data."""

    def __init__(self) -> None:
        """Initialize the proto cluster dictionary."""
        self._clusters = {}

    def __len__(self) -> int:
        """Get the number of proto clusters.
        
        :return: Number of proto clusters.
        :rtype: int
        """        
        return len(self._clusters)

    def add_cluster(self, cluster: ProtoCluster) -> None:
        """Add a proto cluster to the dictionary.
        
        :param cluster: Proto cluster.
        :type cluster: ProtoCluster
        """
        self._clusters[cluster] = cluster

    def assign_gene(self, gene: Gene) -> None:
        """Add a gene to a region.
        
        :param gene: Gene.
        :type gene: Gene
        """
        for cluster_hash in self._clusters.values():
            cluster: ProtoCluster = self._clusters[cluster_hash]
            if cluster.start <= gene.start and gene.end <= cluster.end:
                cluster.genes.append(gene)

    @property
    def clusters(self) -> ty.Dict[int, ProtoCluster]:
        """Get the clusters.
        
        :return: Clusters.
        :rtype: ty.Dict[int, ProtoCluster]
        """
        return self._clusters
    
    def remove_empty_clusters(self) -> None:
        """Remove empty clusters."""
        clusters = self._clusters.copy()
        for cluster in clusters.values():
            if len(cluster) == 0:
                del self._clusters[cluster]


def get_antismash_data(job_id: str) -> ty.Dict[str, ty.Any]:
    """Get AntiSMASH data for a job.
    
    :param job_id: Job ID.
    :type job_id: str
    :return: AntiSMASH data in JSON format.
    :rtype: ty.Dict[str, ty.Any]
    :raises RuntimeError: If less or more than one JSON file is found.
    :raises ValueError: If response status code is not 200.
    """
    logger = logging.getLogger(__name__)

    # Find all JSON files in the directory at job ID url.
    url = os.path.join("https://antismash.secondarymetabolites.org/upload/", job_id)
    response = requests.get(url)

    if response.status_code != 200:
        msg = f"Could not find AntiSMASH data for job ID '{job_id}'."
        logger.error(msg)
        raise ValueError(msg)

    html = response.text
    jsons = []
    for line in html.split("\n"):
        if ".json" in line:
            jsons.append(line.split('"')[1])
    
    # Should have found a single JSON file.
    if len(jsons) != 1:
        msg = "Could not find a single JSON file!"
        logger.error(msg)
        raise RuntimeError(msg)

    # Read contents from that JSON file.
    json_file = jsons[0]
    url = os.path.join("https://antismash.secondarymetabolites.org/upload/", job_id, json_file)
    response = requests.get(url)

    if response.status_code != 200:
        msg = f"Could not find AntiSMASH data for job ID {job_id}!"
        logger.error(msg)
        raise ValueError(msg)

    json_data = response.json()

    return json_data


def get_genes_from_record(record: ty.Dict[str, ty.Any]) -> ty.Dict[str, Gene]:
    """Get genes from an AntiSMASH record.
    
    :param record: AntiSMASH record.
    :type record: ty.Dict[str, ty.Any]
    :return: Genes.
    :rtype: ty.Dict[str, Gene]
    """
    genes = dict()

    features = record["features"]
    for feature in features:
        if feature["type"] != "gene":
            continue
        
        location = feature["location"]
        
        gene_name = feature["qualifiers"].get("gene", None)
        locus_tag = feature["qualifiers"].get("locus_tag", None)

        # Prefer locus tag over gene name.
        if locus_tag is not None:
            name = locus_tag[0]
        else:
            if gene_name is not None:
                name = gene_name[0]
            else:
                raise ValueError("Could not find gene name!")

        if match := re.match(r"\[(\d+):(\d+)\]\(((\+|-))\)", location):
            start = int(match.group(1))
            end = int(match.group(2))
            strand = match.group(3)
        else:
            continue

        gene = Gene(
            name=name,
            start=start,
            end=end,
            strand=strand
        )

        genes[name] = gene
    
    return genes


def get_protoclusters_from_record(
    record: ty.Dict[str, ty.Any],
    categories: ty.List[str] = ["NRPS", "PKS"]
) -> ProtoClusterDict:
    protoclusters = ProtoClusterDict()

    for area in record["areas"]:
        for _, protocluster in area["protoclusters"].items():
            
            category = protocluster["category"]
            product = protocluster["product"]

            if category not in categories:
                continue

            cluster = ProtoCluster(
                start=protocluster["core_start"],
                end=protocluster["core_end"],
                product=product,
                category=category
            )

            protoclusters.add_cluster(cluster)
    
    return protoclusters


def get_modules_from_record(record: ty.Dict[str, ty.Any]) -> ty.Dict[str, ty.List[str]]:
    """Get modules from an AntiSMASH record.

    :param record: AntiSMASH record.
    :type record: ty.Dict[str, ty.Any]
    :return: Modules per gene.
    :rtype: ty.Dict[str, ty.List[str]]
    """
    modules_per_gene = defaultdict(list)  # Modules ordered in the order they appear in the gene.

    modules_type = "antismash.detection.nrps_pks_domains"
    if modules := record["modules"].get(modules_type, None):
        cds_results = modules["cds_results"]
        
        for gene_name, props in cds_results.items():
            predicted_modules = props["modules"]
            if len(predicted_modules) == 0:
                continue

            for predicted_module in predicted_modules:
                module_components = predicted_module["components"]
                for module_component in module_components:
                    hit_id = module_component["domain"]["hit_id"]
                    gene_name = module_component["locus"]
                    modules_per_gene[gene_name].append(hit_id)

    return modules_per_gene


def parse_antismash_json(data: ty.Dict[str, ty.Any]) -> ty.List[ty.Dict[str, ty.Any]]:
    """Parse AntiSMASH JSON data.
    
    :param data: AntiSMASH JSON data.
    :type data: ty.Dict[str, ty.Any]
    :return: Parsed AntiSMASH data.
    :rtype: ty.List[ty.Dict[str, ty.Any]]
    :raises IndexError: If the JSON data is not in the expected format.
    """
    logger = logging.getLogger(__name__)

    # First parse out all genes from the JSON data.
    records = data["records"]

    queries = []
    for record_index, record in enumerate(records):
        dna = record["seq"]["data"]
        genes = get_genes_from_record(record)
        modules = get_modules_from_record(record)
        protoclusters = get_protoclusters_from_record(record)

        # Assign modules to genes.
        for gene_name, gene in genes.items():
            gene.modules = modules.get(gene_name, [])

        # Remove all genes that have no modules.
        genes = {
            name: gene 
            for name, gene 
            in genes.items() 
            if len(gene.modules) > 0
        }

        # Continue if no genes or protoclusters.
        if len(protoclusters) == 0 or len(genes) == 0:
            continue

        # Assign genes to clusters.
        for gene in genes.values():
            protoclusters.assign_gene(gene)

        # Remove clusters without assigned genes.
        protoclusters.remove_empty_clusters()

        logger.debug(f"Found {len(protoclusters)} protoclusters for record {record_index}.")

        for protocluster in protoclusters.clusters.values():
            motif_code_query = protocluster.to_motif_code(dna)
            queries.append(motif_code_query)

    return queries
