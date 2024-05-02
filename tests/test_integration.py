"""Integration tests for the package."""

import json
import os
import typing as ty

import pytest

from retromol.chem import Molecule
from retromol.helpers import TimeoutError, timeout
from retromol.parsing import (
    Result,
    parse_reaction_rules,
    parse_molecular_patterns,
    parse_mol,
)
from retromol.sequencing import parse_modular_natural_product

TESTS_DIR_NAME = os.path.dirname(__file__)
PATH_TO_RULES = os.path.join(TESTS_DIR_NAME, "fixtures", "rules.json")
RULES = json.load(open(PATH_TO_RULES, "r", encoding="utf-8"))
REACTIONS = parse_reaction_rules(json.dumps(RULES["reactions"]))
MONOMERS = parse_molecular_patterns(json.dumps(RULES["monomers"]))


def compile_result(seq: ty.List[str]) -> ty.Dict[str, str]:
    """Compile biosynthetic sequence as strings into a result.

    :param seq: The biosynthetic sequence.
    :type seq: List[str]
    :return: The compiled result.
    :rtype: Dict[str, str]
    """
    items = []

    for motif in seq:
        if motif["identifier"] == "polyketide":
            domains = motif["properties"]["accessory_domains"]
            decoration = motif["properties"]["decoration_type"]

            if set(domains).difference(set(["KR", "DH", "ER"])) == set():
                items.append(f"D{decoration}")
            elif set(domains).difference(set(["KR", "DH"])) == set():
                items.append(f"C{decoration}")
            elif set(domains).difference(set(["KR"])) == set():
                items.append(f"B{decoration}")
            else:
                items.append(f"A{decoration}")

        else:
            raise ValueError(f"Unknown identifier: {motif['identifier']}")

    return items


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
    "test, expected, expected_raises",
    [
        (
            (
                r"CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC"
                r")C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O"
            ),
            ["B7", "B2", "A2", "D7", "B2", "B2"],
            None,
        ),
    ],
)
def test_dummy(test: str, expected: ty.List[str], expected_raises: ty.Any):
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

        sequence = parse_modular_natural_product(result)[0]
        assert compile_result(sequence) == expected
