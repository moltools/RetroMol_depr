
import argparse 
import json
import re
import requests
import time
from collections import defaultdict

from Bio.Seq import Seq

def predict_specificity(src: str):
    url = "https://paras.bioinformatics.nl/api/submit_paras"
    data = {
        "data": {
            "src": src,
            "selectedInputType": "fasta", # or 'genbank'
            "saveActiveSiteSignatures": False,
            "saveExtendedSignatures": False,
            "saveAdenylationDomainSequences": False,
            "selectedSubstrateChoice": "allSubstrates", # or 'commonSubstrates'
            "numPredictionsToReport": 1,
            "useStructureGuidedProfileAlignment": False,
            "firstSeparator": "|",
            "secondSeparator": "_",
            "thirdSeparator": "-"
        }
    }
    headers = {'Content-Type': 'application/json'}
    response = requests.post(url, data=json.dumps(data), headers=headers).json()
    return response

class RegionDictionary:
    def __init__(self):
        self.regions = {}

    def add_region(self, start_region, end_region, category, product):
        self.regions[(start_region, end_region)] = {
            "category": category,
            "product": product,
            "genes": []
        }

    def add_item(self, item, start_item, end_item):
        for region, items in self.regions.items():
            start_region, end_region = region
            if start_region <= start_item and end_item <= end_region:
                items["genes"].append(item)
                break

    def get_regions(self):
        return self.regions
    
    def __len__(self):
        return len(self.regions)

parser = argparse.ArgumentParser()
parser.add_argument("-i", required=True, type=str)
args = parser.parse_args()
path = args.i

# print(path)

with open(args.i, 'r') as f:
    data = json.load(f)

# print(data.keys())
# print(data["version"])
# print(data["input_file"])
# print(data["timings"])
# print(data["taxon"])
# print(data["schema"])
records = data["records"]
# print(len(records))

for record in records:
    # print(record.keys())
    # print(record["id"])
    # print(record["seq"])  # DNA seq. Not needed, very large
    genome = record["seq"]["data"]

    genes = defaultdict(list)
    for feature in record["features"]:
        # print(feature)
        if locus_tag := feature["qualifiers"].get("locus_tag", None):
            if translation := feature["qualifiers"].get("translation", None):
                # print(translation)
                location = feature["location"] # [<start>:<end>](-) or [<start>:<end>](+)
                # print(location, locus_tag, translation)
                if match := re.match(r"\[(\d+):(\d+)\]\(((\+|-))\)", location):
                    start = int(match.group(1))
                    end = int(match.group(2))
                    strand = match.group(3)
                    locus_tag = locus_tag[0]
                    translation = translation[0]
                    genes[locus_tag].append((start, end, strand, translation))

    # print(record["name"])
    # print(record["description"])
    # print(record["dbxrefs"])
    # print(record["annotations"])
    # print(record["letter_annotations"])

    protoclusters = RegionDictionary()  # Probably also want to store what kind of product is predicted to be produced

    for area in record["areas"]:
        # print(area)
        # print(area.keys())  # keys = [start, end, products, protoclusters, candidates, subregions]
        for protocluster_idx, protocluster in area["protoclusters"].items():
            # print(protocluster.keys())  # keys = [category, start, end, core_start, core_end, product, tool]
            if protocluster["category"] == "PKS" or protocluster["category"] == "NRPS":
                # print(protocluster["category"], protocluster["product"], protocluster["start"], protocluster["end"])
                category = protocluster["category"]
                product = protocluster["product"]
                protoclusters.add_region(protocluster["core_start"], protocluster["core_end"], category, product)

    # print(record["modules"].keys())  
    # keys = [
    #   'antismash.detection.hmm_detection', 
    #   'antismash.detection.genefunctions', 
    #   'antismash.detection.nrps_pks_domains',
    #   'antismash.modules.lanthipeptides',
    #   'antismash.modules.lassopeptides',
    #   'antismash.modules.nrps_pks', --> main interest
    #   'antismash.modules.sactipeptides', 
    #   'antismash.modules.t2pks', 
    #   'antismash.modules.tfbs_finder', 
    #   'antismash.modules.thiopeptides', 
    #   'antismash.modules.tta'
    #]
    if modules := record["modules"].get("antismash.detection.nrps_pks_domains", None):
        # print(modules.keys())
        # print(modules["schema_version"])
        # print(modules["record_id"])
        cds_results = modules["cds_results"]
        for locus_tag, props in cds_results.items():
            if gene := genes.get(locus_tag, None):
                for domain in gene:
                    start, end, strand, translation = domain
                    item = {
                        "locus_tag": locus_tag,
                        "start": start,
                        "end": end,
                        "translation": translation,
                        "domains": props
                    }
                    protoclusters.add_item(item, start, end)

        # print(len(protoclusters))
        regions = protoclusters.get_regions()

        for i, (region_loc, region) in enumerate(regions.items()):
            print(f"Region {i} producing {region['product']} ({region['category']}):", region_loc)

            for gene in region["genes"]:
                locus_tag = gene["locus_tag"]

                for domain in genes[locus_tag]:

                    domain_start, domain_end, domain_strand, domain_seq = domain
                    print("\n* " + gene["locus_tag"])
                    for domain in gene["domains"]["domain_hmms"]:
                        print(" >>> ", domain["hit_id"])

                        if locus_tag == "dptA":
                            continue

                        hit_id = domain["hit_id"]
                        if "AMP-binding" in hit_id:
                            dna = Seq(genome[domain_start:domain_end])
                            protein = str(dna.translate())
                            src = f">{region['product']}\n{protein}"
                            response = predict_specificity(src)
                            if payload := response.get("payload", None):
                                if results := payload.get("results", None):
                                    for result in results:
                                        print(">>>>> prediction:", result["data"]["predictions"][0])
                                # time.sleep(1)
    
    # predictions seem to duplicate? every a-domain is predicted 5 times for dptA, dptBC, etc
    # has to do with how the genes are parsed, I think.


