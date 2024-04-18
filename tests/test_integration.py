"""Integration tests for the package."""
import json
import os
import typing as ty

import pytest

from retromol.chem import Molecule
from retromol.helpers import TimeoutError, timeout
from retromol.parsing import Result, parse_reaction_rules, parse_molecular_patterns, parse_mol
from retromol.sequencing import parse_modular_natural_product

PATH_TO_RULES = os.path.join(os.path.dirname(__file__), "fixtures", "rules.json")
RULES = json.load(open(PATH_TO_RULES, "r", encoding="utf-8"))
REACTIONS = parse_reaction_rules(json.dumps(RULES["reactions"]))
MONOMERS = parse_molecular_patterns(json.dumps(RULES["monomers"]))

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
            r"CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O",
            6,
            None
        ),
    ]
)
def test_dummy(test: str, expected: int, expected_raises: ty.Any) -> None:
    """Integration test for parsing and sequencing a natural product compound."""
    molecule = Molecule("test_input", test)
    result = parse_mol_timed(molecule)

    if expected_raises is not None:
        with pytest.raises(expected_raises):
            assert not result.success
    else:
        assert result.success

        sequence = parse_modular_natural_product(result)
        assert len(sequence[0]) == expected
