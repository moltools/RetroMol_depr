[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)
[![Maintainer](https://img.shields.io/badge/Maintainer-davidmeijer-blue)](https://github.com/davidmeijer)
[![Generic badge](https://img.shields.io/badge/Version-0.1.0-green.svg)](https://shields.io/)
![Continuous Integration build in GitHub Actions](https://github.com/moltools/RetroMol/actions/workflows/main.yml/badge.svg?branch=main)


# RetroMol

<img src="./logo.png" alt="logo" width="100">

RetroMol is a Python package for retrosynthetic analysis of molecules.

You can try out RetroMol online [here](https://moltools.bioinformatics.nl/retromol).

## Installation

You can install RetroMol with pip from the root of this repository:
    
```bash
pip install .
```

## Usage CLI

Run the help function for all available commands:

```bash
retromol -h
```

## Run tests 

Run tests:

```bash
python3 -m pytest tests
```