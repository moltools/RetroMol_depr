"""Integration tests for parsing and sequencing modular natural products."""

import json
import os
import typing as ty

import pytest

from retromol.chem import Molecule
from retromol.helpers import timeout
from retromol.parsing import (
    Result,
    parse_reaction_rules,
    parse_molecular_patterns,
    parse_mol,
)

TESTS_DIR_NAME = os.path.dirname(__file__)
PATH_TO_MONOMERS = os.path.join(TESTS_DIR_NAME, "fixtures", "monomers.json")
PATH_TO_REACTIONS = os.path.join(TESTS_DIR_NAME, "fixtures", "reactions.json")
MONOMERS_SRC = json.load(open(PATH_TO_MONOMERS, "r", encoding="utf-8"))
REACTIONS_SRC = json.load(open(PATH_TO_REACTIONS, "r", encoding="utf-8"))
MONOMERS = parse_molecular_patterns(MONOMERS_SRC)
REACTIONS = parse_reaction_rules(REACTIONS_SRC)


@timeout(10)
def parse_mol_timed(mol: Molecule) -> Result:
    """Parse a molecule with a timeout.

    :param mol: The molecule to parse.
    :type mol: Molecule
    :return: The result of the parsing.
    :rtype: Result
    """
    try:
        result = parse_mol(mol, REACTIONS, MONOMERS)
    except Exception:
        result = Result(mol.name, mol.compiled, False)

    return result


@pytest.mark.parametrize(
    "compound_name, test, expected, expected_raises",
    [
        (
            "10-deoxymethynolide",
            r"CC[C@@H]1[C@@H](/C=C/C(=O)[C@@H](C[C@@H]([C@@H]([C@H](C(=O)O1)C)O)C)C)C",
            [
                "polyketide|B2",
                "polyketide|C1",
                "polyketide|A2",
                "polyketide|D2",
                "polyketide|B2",
            ],
            None
        ),
        (
            "13-deoxytedanolide",
            r"C/C=C\[C@H](C)[C@@H]1[C@](O1)(C)[C@H]([C@H]2COC(=O)[C@@H]([C@H]([C@@H](C(=O)[C@@H]([C@H](/C(=C/[C@@H](C(=O)CC[C@H](C2=O)C)C)/C)O)C)C)OC)O)O",
            [
                "polyketide|C2",
                "polyketide|C2",
                "polyketide|B6",
                "polyketide|A2",
                "polyketide|D1",
                "polyketide|A2",
                "polyketide|C2",
                "polyketide|B2",
                "polyketide|A2",
                "polyketide|B5",
            ],
            None
        ),
        (
            "6-deoxyerythronolide_B",
            r"CC[C@@H]1[C@@H]([C@@H]([C@H](C(=O)[C@@H](C[C@@H]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O)C)O)C)C)C)O)C",
            [
                "polyketide|B2",
                "polyketide|B2",
                "polyketide|A2",
                "polyketide|D2",
                "polyketide|B2",
                "polyketide|B2"
            ],
            None
        ),
        (
            "erythromycin",
            r"CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O",
            [
                "polyketide|B7", 
                "polyketide|B2", 
                "polyketide|A2", 
                "polyketide|D7", 
                "polyketide|B2", 
                "polyketide|B2"
            ],
            None
        ),
        (
            "daptomycin",
            r"CCCCCCCCCC(=O)NC(CC1=CNC2=CC=CC=C21)C(=O)NC(CC(=O)N)C(=O)NC(CC(=O)O)C(=O)NC3C(OC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)C(C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C",
            [
                "polyketide|D1",
                "polyketide|D1",
                "polyketide|D1",
                "polyketide|D1",
                "peptide|pubchem|1148",
                "peptide|pubchem|236",
                'peptide|pubchem|424', 
                'peptide|pubchem|205', 
                'peptide|pubchem|750', 
                'peptide|pubchem|389', 
                'peptide|pubchem|424', 
                'peptide|pubchem|602', 
                'peptide|pubchem|424', 
                'peptide|pubchem|750', 
                'peptide|pubchem|617', 
                'peptide|pubchem|237657', 
                'peptide|pubchem|846'
            ],
            None
        )
    ],
)
def test_modular_natural_produts(compound_name: str, test: str, expected: ty.List[str], expected_raises: ty.Any):
    """Integration test for parsing and sequencing a natural product compound.

    :param test: The test input.
    :type test: str
    :param expected: The expected result.
    :type expected: List[str]
    :param expected_raises: The expected exception.
    :type expected_raises: Any
    """
    molecule = Molecule("test_input", test)
    result = parse_mol_timed(molecule)

    if expected_raises is not None:
        with pytest.raises(expected_raises):
            assert not result.success
    else:
        assert result.success

        motif_codes = [seq["motif_code"] for seq in result.sequences]
        assert expected in motif_codes
