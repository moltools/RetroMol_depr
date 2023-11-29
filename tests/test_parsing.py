import json 
import typing as ty
from unittest import TestCase 

from retromol.chem import Molecule 
from retromol.parsing import Result, parse_reaction_rules, parse_molecular_patterns, parse_mol

path_to_rules = "data/rules.json"
rules = json.load(open(path_to_rules, "r"))
reactions = parse_reaction_rules(json.dumps(rules["reactions"]))
monomers = parse_molecular_patterns(json.dumps(rules["monomers"]))

def get_result(mol: Molecule) -> ty.List[ty.List[str]]:
    result = parse_mol(mol, reactions, monomers)
    score = result.score
    return result.biosynthetic_seq, score

class TestParsing10Deoxymethynolide(TestCase):
    def test_parsing_10_deoxymethynolide(self):
        mol = Molecule("10-deoxymethynolide", r"CC[C@@H]1[C@@H](/C=C/C(=O)[C@@H](C[C@@H]([C@@H]([C@H](C(=O)O1)C)O)C)C)C")
        result, score = get_result(mol)
        self.assertTrue(len(result) > 0)
        self.assertTrue(score == 0)
        res_seq = [x[1] for x in result[0]]
        exp_seq = ["B2", "C1", "A2", "D2", "B2"]
        self.assertEqual(res_seq, exp_seq)

class TestParsing6DeoxyerythronolideB(TestCase):
    def test_parsing_6_deoxyerythronoloide_B(self):
        mol = Molecule("6-deoxyerythronolide B", r"CC[C@@H]1[C@@H]([C@@H]([C@H](C(=O)[C@@H](C[C@@H]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O)C)O)C)C)C)O)C")
        result, score = get_result(mol)
        self.assertTrue(len(result) > 0)
        self.assertTrue(score == 0)
        res_seq = [x[1] for x in result[0]]
        exp_seq = ["B2", "B2", "A2", "D2", "B2", "B2"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingAnthracimycin(TestCase):
    def test_parsing_antracimycin(self):
        mol = Molecule("anthracimycin", r"C[C@@H]1/C=C\C=C\[C@H](OC(=O)[C@@H](C(=O)/C=C(/[C@H]2[C@@H]1C=C[C@@H]3[C@@H]2CC=C(C3)C)\O)C)C")
        result, score = get_result(mol)
        self.assertTrue(len(result) > 0)
        self.assertTrue(score == 0)
        res_seq = [x[1] for x in result[0]]
        exp_seq = ["B1", "C1", "C2", "C1", "C1", "D2", "C1", "C1", "A1", "A2"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingBitungolideF(TestCase):
    def test_parsing_bitungolide_F(self):
        mol = Molecule("bitungolide F", r"CC[C@@H]1C=CC(=O)O[C@@H]1[C@H](C)CC[C@H](C[C@@H](/C=C/C=C/C2=CC=CC=C2)O)O")
        result, score = get_result(mol)
        self.assertTrue(len(result) > 0)
        self.assertTrue(score == 0)
        res_seq = [x[1] for x in result[0]]
        exp_seq = ["C1", "C1", "B1", "B1", "D2", "B4", "C1"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingChlorotonilA(TestCase):
    def test_parsing_chlorotonil_A(self):
        mol = Molecule("chlorotonil A", r"C[C@@H]1/C=C\C=C\[C@@H](OC(=O)[C@H](C(=O)C(C(=O)[C@@H]2[C@H]1C=C[C@H]3[C@H]2[C@@H](C=C(C3)C)C)(Cl)Cl)C)C")
        result, score = get_result(mol)
        self.assertTrue(len(result) > 0)
        self.assertTrue(score == 0)
        res_seq = [x[1] for x in result[0]]
        exp_seq = ["B1", "C1", "C2", "C1", "C1", "D2", "C2", "C1", "A1", "A2"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingDaptomycin(TestCase):
    def test_parsing_daptomycin(self):
        mol = Molecule("daptomycin", r"CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)NC3C(OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)C(C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C")
        result, score = get_result(mol)
        self.assertTrue(len(result) > 0)
        self.assertTrue(score == 0)
        res_seq = [x[1] for x in result[0]]
        exp_seq = ["Trp", "Asn", "Asp", "aThr/Thr", "Gly", "Orn", "Asp", "Ala", "Asp", "Gly", "Ser", "3Me-Glu", "Kyn"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingEpothilone(TestCase):
    def test_parsing_epothilone(self):
        mol = Molecule("epothilone", r"C[C@H]1CCC[C@@H]2[C@@H](O2)C[C@H](OC(=O)C[C@H](C(C(=O)[C@@H]([C@H]1O)C)(C)C)O)/C(=C/C3=CSC(=N3)C)/C")
        result, score = get_result(mol)
        self.assertTrue(len(result) > 0)
        self.assertTrue(score == 0)
        res_seq = [x[1] for x in result[0]]
        exp_seq = ["NAc-Cys", "C2", "B1", "C1", "D1", "D2", "B2", "A3", "B1"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingErythromycin(TestCase):
    def test_parsing_erythromycin(self):
        mol = Molecule("erythromycin", r"CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O")
        result, score = get_result(mol)
        self.assertTrue(len(result) > 0)
        self.assertTrue(score == 0)
        res_seq = [x[1] for x in result[0]]
        exp_seq = ["B7", "B2", "A2", "D7", "B2", "B2"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingGephyronicAcid(TestCase):
    def test_parsing_gephyronic_acid(self):
        mol = Molecule("gephyronic acid", r"C[C@@H]1[C@@H](O[C@@](C([C@H]1OC)(C)C)([C@@H](C)C[C@H](C)[C@@H]([C@@]2([C@H](O2)[C@@H](C)C=C(C)C)C)O)O)CC(=O)O")
        result, score = get_result(mol)
        self.assertTrue(len(result) > 0)
        self.assertTrue(score == 0)
        res_seq = [x[1] for x in result[0]]
        exp_seq = ["C2", "C2", "B2", "D2", "B3", "B2", "A1"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingHerboxidiene(TestCase):
    def test_parsing_herboxidiene(self):
        mol = Molecule("herboxidiene", r"C[C@H]1CC[C@@H](O[C@@H]1/C(=C/C=C/[C@@H](C)C[C@@]2([C@H](O2)[C@H](C)[C@H]([C@@H](C)O)OC)C)/C)CC(=O)O")
        result, score = get_result(mol)
        self.assertTrue(len(result) > 0)
        self.assertTrue(score == 0)
        res_seq = [x[1] for x in result[0]]
        exp_seq = ["B2", "C2", "D2", "C1", "C2", "A2", "D1", "C1"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingIndanomycin(TestCase):
    def test_parsing_indanomycin(self):
        mol = Molecule("indanomycin", r"CC[C@H]1CC[C@@H]2[C@@H]1C=C[C@H]([C@H]2C(=O)C3=CC=CN3)/C=C/C=C(\CC)/[C@H]4[C@H](CC[C@@H](O4)[C@@H](C)C(=O)O)C")
        result, score = get_result(mol)
        self.assertTrue(len(result) > 0)
        self.assertTrue(score == 0)
        res_seq = [x[1] for x in result[0]]
        exp_seq = ["A1", "C1", "D4", "C1", "C1", "C1", "C4", "A2", "D1", "C2"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingLactimidomycin(TestCase):
    def test_parsing_lactimidomycin(self):
        mol = Molecule("lactimidomycin", r"C[C@H]1/C=C\C=C\CC/C=C/C(=O)O[C@H]1/C(=C/[C@H](C)C(=O)C[C@@H](CC2CC(=O)NC(=O)C2)O)/C")
        result, score = get_result(mol)
        self.assertTrue(len(result) > 0)
        self.assertTrue(score == 0)
        res_seq = [x[1] for x in result[0]]
        exp_seq = ["Glutarimide", "B1", "A2", "C2", "B2", "C1", "C1", "D1", "C1"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingNarbonolide(TestCase):
    def test_parsing_narbonolide(self):
        mol = Molecule("narbonolide", r"CC[C@@H]1[C@@H](/C=C/C(=O)[C@@H](C[C@@H]([C@@H]([C@H](C(=O)[C@H](C(=O)O1)C)C)O)C)C)C")
        result, score = get_result(mol)
        self.assertTrue(len(result) > 0)
        self.assertTrue(score == 0)
        res_seq = [x[1] for x in result[0]]
        exp_seq = ["B2", "C1", "A2", "D2", "B2", "A2"]
        self.assertEqual(res_seq, exp_seq)