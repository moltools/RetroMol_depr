#!/usr/bin/env python3
"""
Description:    Index parsed ASDB4 proto-clusters.
Usage:          python3 index_asdb.py -i <path to parsed asdb4 JSON dir>
"""
import argparse 
from collections import Counter
import json
from pathlib import Path

COUNTER = Counter()

def cli() -> argparse.Namespace:
    """
    Command line interface.
    
    :return:    Parsed command line arguments.
    :rtype:     argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True, type=str, help="path to parsed asdb4 JSON dir")
    return parser.parse_args()

def parse_file(path: str, num_clusters: int = 0) -> int:
    """
    Parse a single ASDB4 JSON file.
    
    :param str path:    Path to ASDB4 JSON file.
    :return:            Number of clusters in the file.
    :rtype:             int
    """
    data = json.load(open(path))

    # Iterate over each record in the JSON file.
    for record in data:
        proto_clusters = record["proto_clusters"]
        
        # Iterate over each proto-cluster in the record.
        for proto_cluster in proto_clusters:
            num_clusters += 1

            # Iterate over each gene in the proto-cluster.
            for gene in proto_cluster:
                modules = gene["modules"]

                # Every module compiles to a single motif.
                for module in modules:
                    components = module["components"]

                    for component in components:
                        hit_id = component["domain"]["hit_id"]
                        COUNTER[hit_id] += 1

    return num_clusters

def main() -> None:
    """
    Main function.
    """
    # Parse command line arguments.
    args = cli()

    # Get all ASDB4 JSON files.
    file_paths = Path(args.i).rglob("*.json")

    # Iterate over each ASDB4 JSON file.
    num_clusters = 0
    for i, path in enumerate(file_paths):
        try:
            # Parse the file.
            _, num_clusters = parse_file(path, num_clusters)

        except Exception as e:
            # Skip the file if it cannot be parsed.
            print(f"Error parsing {path}: {e}")
            continue
            
        padding = len(str(len(file_paths)))
        print(f"{i}".zfill(padding), end="\r")
    

    print("Keys in proto-cluster JSON files:", COUNTER.keys())
    print("Number of proto-clusters:", num_clusters)
    
    exit(0)

if __name__ == "__main__":
    main()