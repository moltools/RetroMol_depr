#!/usr/bin/env python3
import argparse 
from collections import Counter
import json
from pathlib import Path

"""
num proto clusters = 35014

dict_keys([
    'PKS_KS', 'PKS_AT', 'PP-binding', 'PKS_KR', 'PKS_DH2', 'PKS_DH', 'PKS_ER', 
    'PKS_PP', 'CAL_domain', 'PCP', 'Condensation_LCL', 'AMP-binding', 'Epimerization', 
    'TIGR01720', 'Condensation_DCL', 'Condensation_Dual', 'Thioesterase', 
    'Condensation_Starter', 'TIGR02353', 'Heterocyclization', 'A-OX', 'ACP', 'TD', 
    'cMT', 'Trans-AT_docking', 'PKS_DHt', 'ACP_beta', 'Aminotran_3', 'nMT', 'Cglyc', 
    'X', 'ACPS', 'PT', 'cAT', 'NAD_binding_4', 'Aminotran_1_2', 'FkbH', 'TauD', 'GNAT', 
    'ECH', 'LPG_synthase_C', 'Beta_elim_lyase', 'PS', 'oMT', 'B', 'F'
])
"""

COUNTER = Counter()

class Query:
    pass

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True, type=str, help="path to parsed asdb4 JSON dir")
    return parser.parse_args()

def parse_file(path: str, num_clusters: int = 0) -> dict:
    data = json.load(open(path))

    for record in data:
        proto_clusters = record["proto_clusters"]

        query = []

        for proto_cluster in proto_clusters:
            num_clusters += 1

            for gene in proto_cluster:
                modules = gene["modules"]

                # Every module compiles to a single motif.
                for module in modules:
                    components = module["components"]

                    for component in components:
                        hit_id = component["domain"]["hit_id"]
                        COUNTER[hit_id] += 1

    return {}, num_clusters

def main() -> None:
    args = cli()

    num_clusters = 0
    for i, path in enumerate(Path(args.i).rglob("*.json")):
        try:
            _, num_clusters = parse_file(path, num_clusters)
        except Exception as e:
            print(f"Error parsing {path}: {e}")
            continue
        print(f"{i}".zfill(10), end="\r")
    
    print(COUNTER.keys())
    print(num_clusters)

    exit(0)

if __name__ == "__main__":
    main()