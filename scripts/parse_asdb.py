#!/usr/bin/env python3
"""
Description:    Parse ASDB4 JSON files into a more manageable format.
Usage:          python3 parse_asdb.py -i <path to ASDB4 JSON dir> -o <path to output dir>
"""
import argparse 
import json 
import typing as ty
from pathlib import Path

def cli() -> argparse.Namespace:
    """
    Command line interface.
    
    :return:    Parsed command line arguments.
    :rtype:     argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True, type=str, help="path to JSON dir antiSMASH db.")
    parser.add_argument("-o", required=True, type=str, help="path to output dir.")
    return parser.parse_args()

def parse_gcf(path: str) -> ty.List[ty.Dict[str, ty.Any]]:
    """
    Parse a single ASDB4 JSON file.
    
    :param str path:    Path to ASDB4 JSON file.
    :return:            Proto-clusters in the file.
    :rtype:             ty.List[ty.Dict[str, ty.Any]]
    """
    # Load JSON file.
    data = json.load(open(path))

    # Setup results list.
    results = []

    records = data["records"]
    for record in records:

        # Setup results dict.
        result = dict(
            name = record["name"], 
            description = record["description"], 
            proto_clusters = []
        )

        # Get CDS results.
        cds_results = record["modules"]["antismash.detection.nrps_pks_domains"]["cds_results"]

        # Get region predictions.
        region_predictions = record["modules"]["antismash.modules.nrps_pks"]["region_predictions"]
        for _, regions in region_predictions.items():
            for region in regions:

                ordered_genes = region["ordering"]
                if len(ordered_genes) != 0:
                    
                    proto_cluster = []
                    for gene_accession in ordered_genes:
                        proto_cluster.append(dict(
                            accession = gene_accession,
                            modules = cds_results[gene_accession]["modules"]
                        ))
                    
                    result["proto_clusters"].append(proto_cluster)

        # Only append if there are results:
        if len(result["proto_clusters"]) != 0:
            results.append(result)
    
    return results

def main() -> None:
    args = cli()

    # Create output directory.
    out_dir = Path(args.o)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Get all ASDB4 JSON files.
    file_paths = Path(args.i).rglob("*.json")

    # Iterate over each ASDB4 JSON file.
    failed, total = 0, 0

    # Yield all json files from as path from directory, and parse.
    for i, path in enumerate(file_paths):
        try:
            results = parse_gcf(path)

            # Only write if there are results.
            if len(results) != 0:
                out_path = f"{out_dir}/{path.stem}.json"
                json.dump(results, open(out_path, "w"), indent=4)

                padding = len(str(len(file_paths)))
                print(f'{i}'.zfill(padding), end="\r")

        except:
            failed += 1
        
        total += 1

    # Print failed and total.
    print(f"Failed to parse {failed} out of {total} files.")
    print(f"Output written to {out_dir}.")

    exit(0)

if __name__ == "__main__":
    main()