[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)
[![Maintainer](https://img.shields.io/badge/Maintainer-davidmeijer-blue)](https://github.com/davidmeijer)
[![Generic badge](https://img.shields.io/badge/Version-alpha-green.svg)](https://shields.io/)

# RetroMol

RetroMol is a Python package for retrosynthetic analysis of molecules.

## Installation

Create and activate environment:
    
    ```bash
    conda create -f environment.yml
    conda activate retromol
    ```

## Usage

    ```bash
    retromol -h
    ```

## To do (chronological order)
- [ ] Include stereochemistry (add to monomers and make sure it is recognized, CW, CCW, E, Z, R, S, etc.)
- [ ] Fully parse sample + write tests
- [ ] Include all rules from Alan Healy's rule set
- [ ] Fully pase NPAtlas (PKS/NRPS)
- [ ] Assess performance of rule set
- [ ] Mine for sequence fasta file from NPAtlas results
- [ ] Do some analyses on the fasta sequences
- [ ] Fully parse dataset from Barbara
- [ ] Reintroduce logger