#!/bin/bash

# Check if the required arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_smiles> <output_json> <output_png>"
    exit 1
fi

# Input SMILES and output file
input_smiles="$1"
output_json="$2"
output_png="$3"

# Path to JSON files
reaction_rules_file="./data/reaction_rules.json"
motif_units_file="./data/motif_units.json"
starter_units_file="./data/starter_units.json"
tailoring_units_file="./data/tailoring_units.json"

# Run the retromol command and save output to the specified file
retromol -rr "$reaction_rules_file" \
    -mu "$motif_units_file" \
    -su "$starter_units_file" \
    -tu "$tailoring_units_file" \
    single -i "$input_smiles" > "$output_json" \
    && python3 visualization/visualize_monomer_graph.py -i "$output_json" -o "$output_png"

# Check if the retromol command was successful
if [ $? -ne 0 ]; then
    echo "Error: RetroMol encountered an error."
    exit 1
fi

echo "RetroMol completed successfully."