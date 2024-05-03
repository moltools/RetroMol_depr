#!/usr/bin/env python3
import argparse
import json 
import os

from tqdm import tqdm

from retromol.alignment import (
    ModuleSequence,
    MultipleSequenceAlignment,
    parse_primary_sequence
)
from retromol.parsing import Result


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dir", type=str, required=True, help="Directory containing result JSON files.")
    return parser.parse_args()


def main() -> None:
    args = cli()

    # Get all JSON file paths.
    json_files = [os.path.join(args.dir, x) for x in os.listdir(args.dir) if x.endswith(".json")]

    seqs = []
    # Parse all JSON files.
    for json_path in tqdm(json_files):
        result = Result.from_json(json_path)
        if result.success == True:
            if len(result.sequences) > 0:
                for i, seq in enumerate(result.sequences):
                    i += 1
                    if i > 1:
                        name = f"{result.identifier}_{i}"
                    else:
                        name = result.identifier
                    seq = seq["motif_code"]
                    seq = parse_primary_sequence(seq)
                    seqs.append(ModuleSequence(name, seq))

    # Perform multiple sequence alignment.
    msa = MultipleSequenceAlignment(seqs, 2, 1).get_alignment()

    lens = [len(seq.seq) for seq in msa]
    print(set(lens))

    exit(0)


if __name__ == "__main__":
    main()
