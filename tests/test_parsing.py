from unittest import TestCase 

from retromol.chem import Molecule 
from retromol.parsing import (
    parse_reaction_rules, 
    parse_molecular_patterns, 
    parse_mol
)

RRS = parse_reaction_rules("data/reaction_rules.json")
MUS = parse_molecular_patterns("data/motif_units.json")
SUS = parse_molecular_patterns("data/starter_units.json")
TUS = parse_molecular_patterns("data/tailoring_units.json")

class TestParsingDaptomycin(TestCase):
    def test_parsing_daptomycin(self):
        mol = Molecule("daptomycin", r"CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)NC3C(OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)C(C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C")
        result = parse_mol(mol, RRS, MUS, SUS, TUS)
        res_seq = [x[1] for x in result.biosynthetic_seq[0]]
        exp_seq = ["D1", "D1", "D1", "D1", "Trp", "Asn", "Asp", "aThr/Thr", "Gly", "Orn", "Asp", "Ala", "Asp", "Gly", "Ser", "3Me-Glu", "Kyn"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingErythromycin(TestCase):
    def test_parsing_erythromycin(self):
        mol = Molecule("erythromycin", r"CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O")
        result = parse_mol(mol, RRS, MUS, SUS, TUS)
        res_seq = [x[1] for x in result.biosynthetic_seq[0]]
        exp_seq = ["B7", "B2", "A2", "D7", "B2", "B2"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingGephyronicAcid(TestCase):
    def test_parsing_gephyronic_acid(self):
        mol = Molecule("gephyronic acid", r"C[C@@H]1[C@@H](O[C@@](C([C@H]1OC)(C)C)([C@@H](C)C[C@H](C)[C@@H]([C@@]2([C@H](O2)[C@@H](C)C=C(C)C)C)O)O)CC(=O)O")
        result = parse_mol(mol, RRS, MUS, SUS, TUS)
        res_seq = [x[1] for x in result.biosynthetic_seq[0]]
        exp_seq = ["C2", "C2", "B2", "D2", "B3", "B2", "A1"]
        self.assertEqual(res_seq, exp_seq)
