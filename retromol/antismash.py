
import json
import typing as ty 
import requests


def retrieve_antismash_result(job_id: str, ncbi_acc: str) -> str:
    """Retrieve result JSON from the web.

    :param job_id: The job ID of the AntiSMASH job.
    :type job_id: str
    :param ncbi_acc: The NCBI accession number of the genome.
    :type ncbi_acc: str
    :return: The path to the JSON file.
    :rtype: str
    :raises ValueError: If the request fails.
    """
    url = f"https://antismash.secondarymetabolites.org/upload/{job_id}/{ncbi_acc}.json"

    response = requests.get(url)
    
    if response.status_code == 200:
        return response.text
    
    raise ValueError(f"Failed to retrieve AntiSMASH result: {response.status_code}")


def parse_antismash_result(src: str) -> ty.List[ty.List[str]]:
    """Parse the AntiSMASH result JSON into primary sequences.

    :param src: The path to the JSON file.
    :type src: str
    :return: A list of lists of primary sequences.
    :rtype: List[List[str]]
    """
    data = json.loads(src)

    for record in data["records"]:
        domains = record["modules"]["antismash.detection.nrps_pks_domains"]

        record_id = domains["record_id"]

        for _, domain in domains["cds_results"].items():
            
            # print(domain["modules"])

            modules = domain["modules"]
            print("\n>gene")
            for module in modules:
                components = module["components"]
                print(">domain")
                for component in components:
                    print(component["domain"]["hit_id"])


    # TODO: Parse out modules, and amino acid sequences for A-domains
    # TODO: make predictions with PARAS for A-domains
