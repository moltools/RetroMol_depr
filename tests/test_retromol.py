import json 
import os
import unittest 

from retromol.chem import Molecule 
from retromol.parsing import parse_reaction_rules, parse_molecular_patterns, parse_mol

# Load rules.
path_to_rules_json = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data/rules/rules.json")
rules = json.load(open(path_to_rules_json, "r"))
reactions = parse_reaction_rules(json.dumps(rules["reactions"]))
monomers = parse_molecular_patterns(json.dumps(rules["monomers"]))

# Test parameters.
SCORE_THRESHOLD = 0

def test_score(identifier: str, smiles: str) -> int:
    mol = Molecule(identifier, smiles)
    result = parse_mol(mol, reactions, monomers)
    return result.score

class TestDaptomycin(unittest.TestCase):
    identifier = "daptomycin"
    smiles = r"CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)NC3C(OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)C(C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C"

    def test_score(self) -> None:
        self.assertEqual(test_score(self.identifier, self.smiles), SCORE_THRESHOLD)

class TestErythromycin(unittest.TestCase):
    identifier = "erythromycin"
    smiles = r"CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O"

    def test_score(self) -> None:
        self.assertEqual(test_score(self.identifier, self.smiles), SCORE_THRESHOLD)

if __name__ == "__main__":
    unittest.main()