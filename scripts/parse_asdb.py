#!/usr/bin/env python3
import argparse 
import json 
from pathlib import Path

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True, type=str, help="path to JSON dir antiSMASH db")
    parser.add_argument("-o", required=True, type=str, help="path to output dir")
    return parser.parse_args()

def parse_gcf(path: str) -> dict:
    data = json.load(open(path))

    results = []

    records = data["records"]
    for record in records:

        # Setup results dict.
        result = dict(
            name = record["name"], 
            description = record["description"], 
            proto_clusters = []
        )

        cds_results = record["modules"]["antismash.detection.nrps_pks_domains"]["cds_results"]

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

        results.append(result)
    
    return results

def main() -> None:
    args = cli()

    out_dir = Path(args.o)
    out_dir.mkdir(parents=True, exist_ok=True)

    failed, total = 0, 0

    # Yield all json files from as path from directory, and parse:
    for i, path in enumerate(Path(args.i).rglob("*.json")):
        try:
            results = parse_gcf(path)
            out_path = f"{out_dir}/{path.stem}.json"
            json.dump(results, open(out_path, "w"), indent=4)
            print(f'{i}'.zfill(10), end="\r")

        except:
            failed += 1
        
        total += 1

    print(f"Failed: {failed}/{total}")

    exit(0)

if __name__ == "__main__":
    main()