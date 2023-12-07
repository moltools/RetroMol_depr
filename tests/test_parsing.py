import errno
import json 
import os
import signal
import functools
import typing as ty
from unittest import TestCase 

from retromol.chem import Molecule 
from retromol.helpers import TimeoutError, timeout
from retromol.parsing import parse_reaction_rules, parse_molecular_patterns, parse_mol

path_to_rules = "data/rules.json"
rules = json.load(open(path_to_rules, "r"))
reactions = parse_reaction_rules(json.dumps(rules["reactions"]))
monomers = parse_molecular_patterns(json.dumps(rules["monomers"]))

@timeout(seconds=5)
def get_result(mol: Molecule) -> ty.List[ty.List[str]]:
    result = parse_mol(mol, reactions, monomers)
    score = result.score
    seq = [x[1] for x in result.biosynthetic_seq]
    return seq, score

class TestParsing10Deoxymethynolide(TestCase):
    def test_parsing_10_deoxymethynolide(self):
        mol = Molecule("10-deoxymethynolide", r"CC[C@@H]1[C@@H](/C=C/C(=O)[C@@H](C[C@@H]([C@@H]([C@H](C(=O)O1)C)O)C)C)C")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["B2", "C1", "A2", "D2", "B2"]
        self.assertEqual(res_seq, exp_seq)

class TestParsing13Deoxytedanolide(TestCase):        
    def test_parsing_13_deoxytedanolide(self):
        mol = Molecule("13-deoxytedanolide", r"C/C=C\[C@H](C)[C@@H]1[C@](O1)(C)[C@H]([C@H]2COC(=O)[C@@H]([C@H]([C@@H](C(=O)[C@@H]([C@H](/C(=C/[C@@H](C(=O)CC[C@H](C2=O)C)C)/C)O)C)C)OC)O)O")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["C2", "C2", "B6", "A2", "D1", "A2", "C2", "B2", "A2", "B5"]
        self.assertEqual(res_seq, exp_seq)

class TestParsing6DeoxyerythronolideB(TestCase):
    def test_parsing_6_deoxyerythronoloide_B(self):
        mol = Molecule("6-deoxyerythronolide B", r"CC[C@@H]1[C@@H]([C@@H]([C@H](C(=O)[C@@H](C[C@@H]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O)C)O)C)C)C)O)C")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["B2", "B2", "A2", "D2", "B2", "B2"]
        self.assertEqual(res_seq, exp_seq)

# class TestParsingAbyssomicinC(TestCase):
#     def test_parsing_abyssomicin_C(self):
#         mol = Molecule("abyssomicin C", r"C[C@@H]1C[C@]23OC(=O)C4=C2OC1[C@H](O)C3\C=C\C(=O)[C@@H](C)C[C@@H](C)C4=O")
#         result, score = get_result(mol)
#         self.assertTrue(len(result) > 0)
#         self.assertTrue(score == 0)
#         res_seq = [x[1] for x in result[0]]
#         exp_seq = ["C1", "C1", "C1", "A2", "D2", "A1"]
#         self.assertEqual(res_seq, exp_seq)

class TestParsingAcutiphycin(TestCase):
    def test_parsing_acutiphycin(self):
        mol = Molecule("acutiphycin", r"CCCCC[C@@H]1C/C=C(\[C@H](C(C(=O)[C@H](/C=C(\[C@@H]2C[C@@H](C[C@@](O2)(CC(=O)O1)O)O)/C)C)(C)C)O)/C")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["D1", "D1", "B1", "C2", "B3", "A2", "C2", "A1", "B1", "B1"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingAnthracimycin(TestCase):
    def test_parsing_antracimycin(self):
        mol = Molecule("anthracimycin", r"C[C@@H]1/C=C\C=C\[C@H](OC(=O)[C@@H](C(=O)/C=C(/[C@H]2[C@@H]1C=C[C@@H]3[C@@H]2CC=C(C3)C)\O)C)C")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["B1", "C1", "C2", "C1", "C1", "D2", "C1", "C1", "A1", "A2"]
        self.assertEqual(res_seq, exp_seq)

# class TestParsingAmicoumacin(TestCase):
#     def test_parsing_amicoumacin(self):
#         mol = Molecule("amicoumacin", r"CC(C)C[C@@H](NC(=O)[C@@H](O)[C@@H](O)[C@@H](N)CC(N)=O)C(O)CC(=O)CC(O)CC(=O)CC(=O)O")
#         result, score = get_result(mol)
#         self.assertTrue(len(result) > 0)
#         self.assertTrue(score == 0)
#         res_seq = [x[1] for x in result[0]]
#         exp_seq = ["D", "B5", "Leu", "B1", "A1", "B1", "A1"]
#         self.assertEqual(res_seq, exp_seq)

class TestParsingAmphidinolideJ(TestCase):
    def test_parsing_amphidinolide_J(self):
        mol = Molecule("amphidinolide J", r"CCC/C=C/[C@@H](C)[C@H]1C(/C=C\C([C@H](C=CCCC(=C)[C@H](CC(=O)O1)C)O)C)O")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["D1", "C1", "D9", "B1", "C2", "B1", "C1", "D12", "D8"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingAmphidinolideP(TestCase):
    def test_parsing_amphidinolide_P(self):
        mol = Molecule("amphidinolide P", r"C=C(/C=C/[C@H](O)[C@H](C)C(=C)C)CC=C[C@@H](O)CC(=C)[C@@H](C)C(=O)CC(=O)O")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["D9", "C1", "A8", "C1", "B1", "A9", "A1"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingBitungolideF(TestCase):
    def test_parsing_bitungolide_F(self):
        mol = Molecule("bitungolide F", r"CC[C@@H]1C=CC(=O)O[C@@H]1[C@H](C)CC[C@H](C[C@@H](/C=C/C=C/C2=CC=CC=C2)O)O")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["C1", "C1", "B1", "B1", "D2", "B4", "C1"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingCallystatinA(TestCase):
    def test_parsing_callystatin_A(self):
        mol = Molecule("callystatin A", r"CC[C@H](C)[C@H]([C@H](C)C(=O)[C@H](C)/C=C(\C)/C=C/C[C@@H](C)/C=C(/CC)\C=C\[C@H]1CC=CC(=O)O1)O")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["D2", "B2", "A2", "C2", "C1", "D2", "C4", "C1", "B1", "C1"]
        self.assertEqual(res_seq, exp_seq)

# class TestParsingCarolacton(TestCase):
#     def test_parsing_carolacton(self):
#         mol = Molecule("carolacton", r"C[C@@H]\1CCC[C@@H]([C@H](OC(=O)[C@@H]([C@@H](/C=C1)O)O)/C(=C/[C@@H](C)C(=O)[C@H](C)[C@@H](CC(=O)O)OC)/C)C")
#         result, score = get_result(mol)
#         self.assertTrue(len(result) > 0)
#         self.assertTrue(score == 0)
#         res_seq = [x[1] for x in result[0]]
#         exp_seq = [] # Has multiple output sequences, probably...
#         self.assertEqual(res_seq, exp_seq)

class TestParsingCurvularideC(TestCase):
    def test_parsing_curvularide_C(self):
        mol = Molecule("curvularide C", r"CC[C@H](C)[C@@H](CO)NC(=O)/C=C/[C@](C)([C@H]([C@@H](C)C[C@@H](CC)O)O)OC")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["D5", "D2", "B7", "C1", "Ile/aIle"]
        self.assertEqual(res_seq, exp_seq)

# class TestParsingChaetoglobosinA(TestCase):
#     def test_parsing_chaetoglobosin_A(self):
#         mol = Molecule("chaetoglobosin A", r"C[C@H]\1C/C=C/[C@H]2[C@H]3[C@](O3)([C@H]([C@@H]4[C@@]2(C(=O)/C=C/C(=O)[C@@H](/C(=C1)/C)O)C(=O)N[C@H]4CC5=CNC6=CC=CC=C65)C)C")
#         result, score = get_result(mol)
#         self.assertTrue(len(result) > 0)
#         self.assertTrue(score == 0)
#         res_seq = [x[1] for x in result[0]]
#         exp_seq = ["C2", "C1", "C1", "D2", "C2", "B11", "C1", "A", "Trp"]
#         self.assertEqual(res_seq, exp_seq)

class TestParsingChlorotonilA(TestCase):
    def test_parsing_chlorotonil_A(self):
        mol = Molecule("chlorotonil A", r"C[C@@H]1/C=C\C=C\[C@@H](OC(=O)[C@H](C(=O)C(C(=O)[C@@H]2[C@H]1C=C[C@H]3[C@H]2[C@@H](C=C(C3)C)C)(Cl)Cl)C)C")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["B1", "C1", "C2", "C1", "C1", "D2", "C2", "C1", "A1", "A2"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingDaptomycin(TestCase):
    def test_parsing_daptomycin(self):
        mol = Molecule("daptomycin", r"CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)NC3C(OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)C(C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["Trp", "Asn", "Asp", "aThr/Thr", "Gly", "Orn", "Asp", "Ala", "Asp", "Gly", "Ser", "3Me-Glu", "Kyn"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingDictyostatin(TestCase):
    def test_parsing_dictyostatin(self):
        mol = Molecule("dictyostatin", r"C[C@H]1CC[C@H]([C@@H]([C@@H](OC(=O)/C=C\C=C\[C@H]([C@H](C[C@@H](/C=C\[C@@H]([C@@H]([C@H](C1)C)O)C)O)O)C)[C@@H](C)/C=C\C=C)C)O")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["C1", "C2", "B2", "B1", "D2", "D2", "B2", "C1", "B1", "B2", "C1", "C1"]
        self.assertEqual(res_seq, exp_seq)

# class TestParsingDiscodermolide(TestCase):
#     def test_parsing_discodermolide(self):
#         mol = Molecule("discodermolide", r"C[C@H]1[C@@H](OC(=O)[C@@H]([C@H]1O)C)C[C@@H](/C=C\[C@H](C)[C@@H]([C@@H](C)/C=C(/C)\C[C@H](C)[C@H]([C@H](C)[C@H]([C@@H](C)/C=C\C=C)OC(=O)N)O)O)O")
#         result, score = get_result(mol)
#         self.assertTrue(len(result) > 0)
#         self.assertTrue(score == 0)
#         res_seq = [x[1] for x in result[0]]
#         exp_seq = ["C1", "C2", "B2", "B2", "D2", "C2", "B2", "C1", "B1", "B2", "B2"]
#         self.assertEqual(res_seq, exp_seq)

class TestParsingEpothilone(TestCase):
    def test_parsing_epothilone(self):
        mol = Molecule("epothilone", r"C[C@H]1CCC[C@@H]2[C@@H](O2)C[C@H](OC(=O)C[C@H](C(C(=O)[C@@H]([C@H]1O)C)(C)C)O)/C(=C/C3=CSC(=N3)C)/C")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["NAc-Cys", "C2", "B1", "C1", "D1", "D2", "B2", "A3", "B1"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingErythromycin(TestCase):
    def test_parsing_erythromycin(self):
        mol = Molecule("erythromycin", r"CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["B7", "B2", "A2", "D7", "B2", "B2"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingGephyronicAcid(TestCase):
    def test_parsing_gephyronic_acid(self):
        mol = Molecule("gephyronic acid", r"C[C@@H]1[C@@H](O[C@@](C([C@H]1OC)(C)C)([C@@H](C)C[C@H](C)[C@@H]([C@@]2([C@H](O2)[C@@H](C)C=C(C)C)C)O)O)CC(=O)O")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["C2", "C2", "B2", "D2", "B3", "B2", "A1"]
        self.assertEqual(res_seq, exp_seq)

# class TestParsingHarzianicAcid(TestCase):
#     def test_parsing_harzianic_acid(self):
#         mol = Molecule("harzianic acid", r"CCC/C=C/C=C/C(=C\1/C(=O)C(N(C1=O)C)CC(C(C)C)(C(=O)O)O)/O")
#         result, score = get_result(mol)
#         self.assertTrue(len(result) > 0)
#         self.assertTrue(score == 0)
#         res_seq = [x[1] for x in result[0]]
#         exp_seq = ["D1", "C1", "C1", "A1", "AA"]
#         self.assertEqual(res_seq, exp_seq)

class TestParsingHerboxidiene(TestCase):
    def test_parsing_herboxidiene(self):
        mol = Molecule("herboxidiene", r"C[C@H]1CC[C@@H](O[C@@H]1/C(=C/C=C/[C@@H](C)C[C@@]2([C@H](O2)[C@H](C)[C@H]([C@@H](C)O)OC)C)/C)CC(=O)O")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["B2", "C2", "D2", "C1", "C2", "A2", "D1", "C1"]
        self.assertEqual(res_seq, exp_seq)

# class TestParsingHymenosetin(TestCase):
#     def test_parsing_hymenosetin(self):
#         mol = Molecule("hymenosetin", r"C/C=C/[C@@H]1C(=C[C@@H]2C[C@@H](CC[C@H]2[C@]1(C)/C(=C\3/C(=O)[C@@H](NC3=O)[C@@H](C)O)/O)C)C")
#         result, score = get_result(mol)
#         self.assertTrue(len(result) > 0)
#         self.assertTrue(score == 0)
#         res_seq = [x[1] for x in result[0]]
#         exp_seq = ["C1", "C2", "C1", "D2", "D1", "C2", "A1", "aThr/Thr"]
#         self.assertEqual(res_seq, exp_seq)

class TestParsingIndanomycin(TestCase):
    def test_parsing_indanomycin(self):
        mol = Molecule("indanomycin", r"CC[C@H]1CC[C@@H]2[C@@H]1C=C[C@H]([C@H]2C(=O)C3=CC=CN3)/C=C/C=C(\CC)/[C@H]4[C@H](CC[C@@H](O4)[C@@H](C)C(=O)O)C")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["A1", "C1", "D4", "C1", "C1", "C1", "C4", "A2", "D1", "C2"]
        self.assertEqual(res_seq, exp_seq)

# class TestParsingIrcinianin(TestCase):
#     def test_parsing_ircinianin(self):
#         mol = Molecule("ircinianin", r"CC(C=C/C=C(\C)CCCc1ccoc1)=CCC[C@H](C)C=C1OC(=O)C(C)=C1O")
#         result, score = get_result(mol)
#         self.assertTrue(len(result) > 0)
#         self.assertTrue(score == 0)
#         res_seq = [x[1] for x in result[0]]
#         exp_seq = ["D1", "D2", "C1", "C2", "C1", "D2", "D11", "A2"]
#         self.assertEqual(res_seq, exp_seq)

# class TestParsingIriomoteolide1a(TestCase):
#     def test_parsing_iriomoteolide1a(self):
#         mol = Molecule("iriomoteolide 1a", r"C[C@H]1C/C=C/[C@@]([C@@]2(CC(=C)C[C@@H](O2)C/C=C/[C@@H]([C@@H](/C(=C\C(=O)O[C@@H]1C[C@H](C)[C@H](C)O)/C)C)O)O)(C)O")
#         result, score = get_result(mol)
#         self.assertTrue(len(result) > 0)
#         self.assertTrue(score == 0)
#         res_seq = [x[1] for x in result[0]]
#         exp_seq = ["D8", "B2", "D1", "C", "A1", "A8", "B1", "C1", "B2", "C6"]
#         self.assertEqual(res_seq, exp_seq)

# class TestParsingIriomoteolide3a(TestCase):
#     def test_parsing_iriomoteolide3a(self):
#         mol = Molecule("iriomoteolide 1a", r"C/C=C/C/C=C/[C@H](C)C[C@@H]([C@@H]1C[C@H]2[C@@H](O2)/C=C/[C@@H]([C@H](/C=C/C[C@@H](CC(=O)O1)C)O)O)O")
#         result, score = get_result(mol)
#         self.assertTrue(len(result) > 0)
#         self.assertTrue(score == 0)
#         res_seq = [x[1] for x in result[0]]
#         exp_seq = ["C1", "C1", "D8", "B5", "D1", "C1", "C", "B1", "C1", "D8"]
#         self.assertEqual(res_seq, exp_seq)

# class TestParsingJerangolidA(TestCase):
#     def test_linearization_jerangolid_A(self):
#         mol = Molecule(r"CC[C@@H]1C(=CC[C@@H](O1)/C(=C/[C@H](C)/C=C/[C@H]2CC(=C(C(=O)O2)CO)OC)/C)C")
#         results = linearize(mol)
#         expected = Molecule(r"CC=CC(C)=CC[C@@H](O)/C(C)=C/[C@H](C)/C=C/[C@H](O)CC(O)=C(CO)C(=O)O")
#         self.assertEqual(len(results), 1)
#         self.assertTrue(results[0] == expected)

#     def test_sequencing_jerangolid_A(self):
#         mol = Molecule(r"CC=CC(C)=CC[C@@H](O)/C(C)=C/[C@H](C)/C=C/[C@H](O)CC(O)=C(CO)C(=O)O")
#         result = [str(u) for u in sequence(mol)]
#         expected = ["C2", "C1", "B2", "C2", "C1", "B1", "A6"]
#         self.assertEqual(result, expected)

# class TestParsingKirromycin(TestCase):
#     def test_linearization_kirromycin(self):
#         mol = Molecule(r"CC[C@H](C(=O)NC/C=C/C=C(\C)/[C@H]([C@@H](C)[C@H]1[C@H]([C@H]([C@H](O1)/C=C/C=C/C=C(\C)/C(=O)C2=C(C=CNC2=O)O)O)O)OC)[C@@]3([C@@H]([C@@H](C([C@@H](O3)/C=C/C=C/C)(C)C)O)O)O")
#         results = linearize(mol)
#         expected = Molecule(r"C/C=C/C=C/[C@H](O)C(C)(C)[C@@H](O)[C@@H](O)C(=O)[C@H](CC)C(=O)NC/C=C/C=C(\C)[C@@H](O)C(C)=C[C@@H](O)[C@@H](O)[C@H](O)/C=C/C=C/C=C(\C)C(=O)CC(=O)NCC(=O)O")
#         self.assertEqual(len(results), 1)
#         self.assertTrue(results[0] == expected)

#     def test_sequencing_kirromycin(self):
#         mol = Molecule(r"C/C=C/C=C/[C@H](O)C(C)(C)[C@@H](O)[C@@H](O)C(=O)[C@H](CC)C(=O)NC/C=C/C=C(\C)[C@@H](O)C(C)=C[C@@H](O)[C@@H](O)[C@H](O)/C=C/C=C/C=C(\C)C(=O)CC(=O)NCC(=O)O")
#         result = [str(u) for u in sequence(mol)]
#         expected = ["C1", "C1", "B3", "B5", "A4", "Gly", "C1", "C2", "B2", "C", "B5", "C1", "C1", "C2", "A1", "Gly"] # TODO: check correctness test.
#         self.assertEqual(result, expected)

class TestParsingLactimidomycin(TestCase):
    def test_parsing_lactimidomycin(self):
        mol = Molecule("lactimidomycin", r"C[C@H]1/C=C\C=C\CC/C=C/C(=O)O[C@H]1/C(=C/[C@H](C)C(=O)C[C@@H](CC2CC(=O)NC(=O)C2)O)/C")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["Glutarimide", "B1", "A2", "C2", "B2", "C1", "C1", "D1", "C1"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingLankamycin(TestCase):
    def test_parsing_lankamycin(self):
        mol = Molecule("lankamycin", r"C[C@H]([C@@H](O)[C@@H](C)C[C@](C)(O)C(=O)[C@H](C)[C@@H](O)[C@@H](C)[C@H](O)[C@@H](C)[C@H](C)O)[C@H](O)[C@@H](C)C(=O)O")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["B2", "B2", "B2", "A7", "D2", "B2", "B2"]
        self.assertEqual(res_seq, exp_seq)

# class TestParsingLatrunculinA(TestCase):
#     def test_linearization_latrunculin_A(self):
#         mol = Molecule(r"C[C@H]/1CC[C@@H]2C[C@H](C[C@@](O2)([C@@H]3CSC(=O)N3)O)OC(=O)/C=C(\CC/C=C/C=C1)/C")
#         results = linearize(mol)
#         expected = Molecule(r"C/C(=C/C(=O)O)CC/C=C/C=C\[C@@H](C)CC[C@@H](O)C[C@@H](O)CC(=O)[C@@H](N)CS")
#         self.assertEqual(len(results), 1)
#         self.assertTrue(results[0] == expected)

#     def test_sequencing_latrunculin_A(self):
#         mol = Molecule(r"C/C(=C/C(=O)O)CC/C=C/C=C\[C@@H](C)CC[C@@H](O)C[C@@H](O)CC(=O)[C@@H](N)CS")
#         result = [str(u) for u in sequence(mol)]
#         expected = ["NRP", "A1", "B1", "B1", "D2", "C1", "C1", "D1", "C6"]
#         self.assertEqual(result, expected)

# class TestParsingLeiodermatolide(TestCase):
#     def test_linearization_leiodermatolide(self):
#         mol = Molecule(r"CC[C@H]1[C@@H]([C@@](CC(=O)O1)(C/C=C/C=C(\C)/[C@@H]2[C@H](/C=C\C=C/[C@@H]([C@H]([C@H]([C@@H](/C(=C/CCC(=O)O2)/C)C)O)C)OC(=O)N)C)O)C")
#         results = linearize(mol)
#         expected = Molecule(r"CC[C@H](O)[C@H](C)[C@](O)(C/C=C/C=C(\C)[C@@H](O)[C@@H](C)/C=C\C=C/[C@H](O)[C@@H](C)[C@@H](O)[C@H](C)/C(C)=C/CCC(=O)O)CC(=O)O")
#         self.assertEqual(len(results), 1)
#         self.assertTrue(results[0] == expected)

#     def test_sequencing_leiodermatolide(self):
#         mol = Molecule(r"CC[C@H](O)[C@H](C)[C@](O)(C/C=C/C=C(\C)[C@@H](O)[C@@H](C)/C=C\C=C/[C@H](O)[C@@H](C)[C@@H](O)[C@H](C)/C(C)=C/CCC(=O)O)CC(=O)O")
#         self.assertRaises(SequencingError, sequence, mol) # Has two starting points (R-COOH) groups.

# class TestParsingLovastatin(TestCase):
#     def test_linearization_lovastatin(self):
#         mol = Molecule(r"CC[C@H](C)C(=O)O[C@H]1C[C@H](C=C2[C@H]1[C@H]([C@H](C=C2)C)CC[C@@H]3C[C@H](CC(=O)O3)O)C")
#         results = linearize(mol)
#         expected = Molecule(r"")
#         self.assertEqual(len(results), 1)
#         self.assertTrue(results[0] == expected)

#     def test_sequencing_lovastatin(self):
#         mol = Molecule(r"")
#         result = [str(u) for u in sequence(mol)]
#         expected = []
#         self.assertEqual(result, expected)

class TestParsingMacrolactin(TestCase):
    def test_parsing_macrolactin(self):
        mol = Molecule("macrolactin", r"CC1CCC/C=C/C=C/C(CC(C/C=C\C=C\C(C/C=C/C=C\C(=O)O1)O)O)O")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["B1", "D1", "C1", "C1", "B1", "B1", "C1", "C1", "B1", "C1", "C1"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingMegalomycinA(TestCase):
    def test_parsing_megalomycin_A(self):
        mol = Molecule("megalomycin A", r"CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)O)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O[C@H]4C[C@H]([C@H]([C@@H](O4)C)O)N(C)C)C)C)O)(C)O")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["B7", "B2", "A2", "D7", "B2", "B2"]
        self.assertEqual(res_seq, exp_seq)

# class TestParsingMicacocidinA(TestCase):
#     def test_linearization_micacocidin_A(self):
#         mol = Molecule(r"CCCCCC1=C(C(=CC=C1)[O-])C2=N[C@H](CS2)[C@@H]3N([C@@H](CS3)[C@@H](C(C)(C)C4=N[C@@](CS4)(C)C(=O)[O-])O)C")
#         results = linearize(mol)
#         expected = Molecule(r"CCCCCc1cccc([O-])c1C(=O)N[C@H](CS)C(=O)N[C@@H](CS)[C@H](O)C(C)(C)C(=O)N[C@](C)(CS)C(=O)[O-]")
#         self.assertEqual(len(results), 1)
#         self.assertTrue(results[0] == expected)
    
#     def test_sequencing_micacocidin_A(self):
#         mol = Molecule(r"CCCCCc1cccc([O-])c1C(=O)N[C@H](CS)C(=O)N[C@@H](CS)[C@H](O)C(C)(C)C(=O)N[C@](C)(CS)C(=O)[O-]")
#         result = [str(u) for u in sequence(mol)]
#         expected = ["Cys", "Cys", "B3", "AlphaAminoAcid"]
#         self.assertEqual(result, expected)

class TestParsingMigrastatin(TestCase):
    def test_parsing_migrastatin(self):
        mol = Molecule("migrastatin", r"C[C@@H]1/C=C(\[C@H](OC(=O)/C=C/CC/C=C/[C@@H]([C@H]1O)OC)[C@H](C)C(=O)CCCC2CC(=O)NC(=O)C2)/C")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["Glutarimide", "D1", "A2", "B2", "C2", "B5", "C1", "D1", "C1"]
        self.assertEqual(res_seq, exp_seq)

class TestParsingNarbonolide(TestCase):
    def test_parsing_narbonolide(self):
        mol = Molecule("narbonolide", r"CC[C@@H]1[C@@H](/C=C/C(=O)[C@@H](C[C@@H]([C@@H]([C@H](C(=O)[C@H](C(=O)O1)C)C)O)C)C)C")
        res_seq, score = get_result(mol)
        self.assertTrue(score == 0)
        exp_seq = ["B2", "C1", "A2", "D2", "B2", "A2"]
        self.assertEqual(res_seq, exp_seq)

# class TestParsingPederin(TestCase):
#     def test_linearization_pederin(self):
#         mol = Molecule(r"C[C@H]1[C@H](O[C@](CC1=C)([C@@H](C(=O)N[C@H]([C@@H]2C[C@H](C([C@H](O2)C[C@@H](COC)OC)(C)C)O)OC)O)OC)C")
#         results = linearize(mol)
#         # Undirected ether linearization can lead to two different results.
#         expected1 = Molecule(r"C=C(CC(=O)[C@H](O)C(=O)N[C@@H](O)[C@@H](O)C[C@@H](O)C(C)(C)C=C[C@H](O)CO)[C@@H](C)[C@@H](C)O")
#         expected2 = Molecule(r"C=C(CC(=O)[C@H](O)C(=O)NC(O)=CC[C@@H](O)C(C)(C)[C@H](O)C[C@H](O)CO)[C@@H](C)[C@@H](C)O")
#         self.assertEqual(len(results), 2)
#         self.assertTrue(results[0] == expected1 or results[0] == expected2)
#         self.assertTrue(results[1] == expected1 or results[1] == expected2)

#     def test_sequencing_pederin(self):
#         mol1 = Molecule(r"C=C(CC(=O)[C@H](O)C(=O)N[C@@H](O)[C@@H](O)C[C@@H](O)C(C)(C)C=C[C@H](O)CO)[C@@H](C)[C@@H](C)O")
#         mol2 = Molecule(r"C=C(CC(=O)[C@H](O)C(=O)NC(O)=CC[C@@H](O)C(C)(C)[C@H](O)C[C@H](O)CO)[C@@H](C)[C@@H](C)O")
#         result1 = [str(u) for u in sequence(mol1)]
#         result2 = [str(u) for u in sequence(mol2)]
#         expected1 = ["B2", "A8", "A5", "AlphaAminoAcid", "B1", "B3", "C1"]
#         expected2 = ["B2", "A8", "A5", "AlphaAminoAcid", "C1", "B3", "B1"]
#         self.assertEqual(result1, expected1)
#         self.assertEqual(result2, expected2)

# class TestParsingPelorusideA(TestCase):
#     def test_linearization_peloruside_A(self):
#         mol = Molecule(r"CC[C@@H](CO)/C=C(/C)\[C@@H]1C[C@H](C[C@@H](C([C@@]2([C@@H]([C@@H](C[C@@H](O2)C[C@H]([C@@H](C(=O)O1)O)OC)OC)O)O)(C)C)O)OC")
#         results = linearize(mol)
#         expected = Molecule(r"CC[C@H](/C=C(/C)[C@@H](O)C[C@@H](O)C[C@H](O)C(C)(C)C(=O)[C@H](O)[C@H](O)C[C@@H](O)C[C@@H](O)[C@H](O)C(=O)O)CO")
#         self.assertEqual(len(results), 1)
#         self.assertTrue(results[0] == expected)

#     def test_sequencing_peloruside_A(self):
#         mol = Molecule(r"CC[C@H](/C=C(/C)[C@@H](O)C[C@@H](O)C[C@H](O)C(C)(C)C(=O)[C@H](O)[C@H](O)C[C@@H](O)C[C@@H](O)[C@H](O)C(=O)O)CO")
#         result = [str(u) for u in sequence(mol)]
#         expected = ["D6", "C2", "B1", "B1", "B3", "A5", "B1", "B1", "B5"]
#         self.assertEqual(result, expected)

# class TestParsingPenicillinG(TestCase):
#     def test_parsing_penicillin_G(self):
#         mol = Molecule("penicillin G", r"CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C")
#         result, score = get_result(mol)
#         self.assertTrue(len(result) > 0)
#         self.assertTrue(score == 0)
#         res_seq = [x[1] for x in result[0]]
#         exp_seq = []
#         self.assertEqual(res_seq, exp_seq)

# # class TestParsingPericoniasinI(TestCase):
# #     def test_linearization_periconiasin_I(self):
# #         mol = Molecule(r"C/C/1=C/C[C@@H](CC(=O)[C@]23[C@@H](C1)[C@H](C(=C([C@H]2[C@@H](NC3=O)CC(C)C)C)C)O)O")
# #         results = linearize(mol)
# #         expected = Molecule(r"")
# #         self.assertEqual(len(results), 1)
# #         self.assertTrue(results[0] == expected)

# # def test_sequencing_periconiasin_I(self):
# #     mol = Molecule(r"")
# #     result = [str(u) for u in sequence(mol)]
# #     expected = []
# #     self.assertEqual(result, expected)

# class TestParsingRatjadon(TestCase):
#     def test_linearization_ratjadon(self):
#         mol = Molecule(r"C/C=C/[C@H]1[C@H]([C@@H](C[C@H](O1)[C@@H](/C=C/C=C(\C)/C[C@@H](C)/C=C(/C)\C=C\[C@H]2CC=CC(=O)O2)O)O)C")
#         results = linearize(mol)
#         expected = Molecule(r"C/C=C/[C@H](O)[C@@H](C)[C@H](O)CC=C(O)/C=C/C=C(\C)C[C@@H](C)/C=C(C)\C=C\[C@H](O)CC=CC(=O)O")
#         self.assertEqual(len(results), 1)
#         self.assertTrue(results[0] == expected)

#     def test_sequencing_ratjadon(self):
#         mol = Molecule(r"C/C=C/[C@H](O)[C@@H](C)[C@H](O)CC=C(O)/C=C/C=C(\C)C[C@@H](C)/C=C(C)\C=C\[C@H](O)CC=CC(=O)O")
#         result = [str(u) for u in sequence(mol)]
#         expected = ["C1", "B2", "B1", "D11", "C1", "C2", "D2", "C2", "C1", "B1", "C1"]
#         self.assertEqual(result, expected)

# class TestParsingSalinosporimide(TestCase):
#     def test_linearization_salinosporimide(self):
#         mol = Molecule(r"CCCC1C(=O)NC2(C1(OC2=O)C)C(C3CCCC=C3)O")
#         results = linearize(mol)   
#         expected = Molecule(r"CCCCC(=O)NC(C(=O)OCC)C(O)C1C=CCCC1")
#         self.assertEqual(len(results), 1)
#         self.assertTrue(results[0] == expected)

#     # def test_sequencing_salinosporimide(self):
#     #     mol = Molecule(r"CCCCC(=O)NC(C(=O)OCC)C(O)C1C=CCCC1")
#     #     result = [str(u) for u in sequence(mol)]
#     #     expected = []
#     #     self.assertEqual(result, expected)

# class TestParsingSoraphenA(TestCase):
#     def test_linearization_soraphen_A(self):
#         mol = Molecule(r"C[C@H]1/C=C/[C@H]([C@H](CCCC[C@H](OC(=O)[C@H]([C@@]2([C@@H]([C@H]([C@@H]([C@H]1O2)C)O)OC)O)C)C3=CC=CC=C3)OC)OC")
#         results = linearize(mol)
#         expected = Molecule(r"C[C@H](C(=O)O)C(=O)[C@H](O)[C@@H](O)[C@H](C)[C@@H](O)[C@@H](C)/C=C/[C@@H](O)[C@@H](O)CCCC[C@H](O)c1ccccc1")
#         self.assertEqual(len(results), 1)
#         self.assertTrue(results[0] == expected)

#     def test_sequencing_soraphen_A(self):
#         mol = Molecule(r"C[C@H](C(=O)O)C(=O)[C@H](O)[C@@H](O)[C@H](C)[C@@H](O)[C@@H](C)/C=C/[C@@H](O)[C@@H](O)CCCC[C@H](O)c1ccccc1")
#         result = [str(u) for u in sequence(mol)]
#         expected = ["B1", "D1", "D5", "B1", "C2", "B2", "B5", "A2"]
#         self.assertEqual(result, expected)

# class TestParsingSpiculoicAcidA(TestCase):
#     def test_linearization_spiculoic_acid_A(self):
#         mol = Molecule(r"CC[C@@H]1[C@@H]2[C@@H]([C@H](C1=O)C)C(=C[C@@]([C@@]2(CC)C(=O)O)(CC)/C=C/C3=CC=CC=C3)CC")
#         results = linearize(mol)
#         expected = Molecule(r"CCC(=C[C@@H](C)C(=O)[C@@H](C=C(CC)C(=O)O)CC)C=C(/C=C/c1ccccc1)CC")
#         self.assertEqual(len(results), 1)
#         self.assertTrue(results[0] == expected)

#     def test_sequencing_spiculoic_acid_A(self):
#         mol = Molecule(r"CCC(=C[C@@H](C)C(=O)[C@@H](C=C(CC)C(=O)O)CC)C=C(/C=C/c1ccccc1)CC")
#         result = [str(u) for u in sequence(mol)]
#         expected = ["C4", "C4", "C2", "A4", "C4"]
#         self.assertEqual(result, expected)

# class TestParsingSpongidepsin(TestCase):
#     def test_linearization_spongidepsin(self):
#         mol = Molecule(r"C[C@H]1CCC(CC(OC(=O)[C@@H](N(C(=O)[C@H](C1)C)C)CC2=CC=CC=C2)CCCC#C)C")
#         results = linearize(mol)
#         expected = Molecule(r"C#CCCCC(O)CC(C)CC[C@H](C)C[C@H](C)C(=O)N[C@@H](Cc1ccccc1)C(=O)O")
#         self.assertEqual(len(results), 1)
#         self.assertTrue(results[0] == expected)

#     def test_sequencing_spongidepsin(self):
#         mol = Molecule(r"C#CCCCC(O)CC(C)CC[C@H](C)C[C@H](C)C(=O)N[C@@H](Cc1ccccc1)C(=O)O")
#         result = [str(u) for u in sequence(mol)]
#         expected = ["D1", "B1", "D8", "D2", "D2", "Phe"]
#         self.assertEqual(result, expected)

# class TestParsingThailanstatinA(TestCase):
#     def test_linearization_thailanstatin(self):
#         mol = Molecule(r"C[C@H]1C[C@H]([C@H](O[C@H]1C/C=C(\C)/C=C/[C@@H]2[C@H]([C@@]3(C[C@H](O2)CC(=O)O)CO3)O)C)NC(=O)/C=C\[C@H](C)OC(=O)C")
#         results = linearize(mol)
#         expected = Molecule(r"C=C(CC=CC(=O)O)[C@H](O)[C@H](O)/C=C/C(C)=C/C=C[C@@H](C)C[C@@H](NC(=O)/C=C\[C@H](C)O)[C@@H](C)O")
#         self.assertEqual(len(results), 1)
#         self.assertTrue(results[0] == expected)

#     def test_sequencing_thailanstatin(self):
#         mol = Molecule(r"C=C(CC=CC(=O)O)[C@H](O)[C@H](O)/C=C/C(C)=C/C=C[C@@H](C)C[C@@H](NC(=O)/C=C\[C@H](C)O)[C@@H](C)O")
#         result = [str(u) for u in sequence(mol)]
#         expected = ["C1", "aThr/Thr", "D2", "C1", "C2", "C1", "B5", "A8", "C1"]
#         self.assertEqual(result, expected)
    
# class TestParsingTheopederinA(TestCase):
#     def test_linearization_theopederin_A(self):
#         mol = Molecule(r"C=C(CC(=O)[C@@H](O)C(=O)N[C@H]1OCO[C@H]([C@@H](O)C(C)(C)C=C[C@H]2CCCC(O)O2)[C@@H]1O)[C@@H](C)[C@@H](C)O")
#         results = linearize(mol)
#         expected = Molecule(r"C=C(CC(=O)[C@@H](O)C(=O)N[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(C)(C)C=C[C@H](O)CCCC=O)[C@@H](C)[C@@H](C)O")
#         self.assertEqual(len(results), 1)
#         self.assertTrue(results[0] == expected)

#     def test_sequencing_theopederin_A(self):
#         mol = Molecule(r"C=C(CC(=O)[C@@H](O)C(=O)N[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(C)(C)C=C[C@H](O)CCCC=O)[C@@H](C)[C@@H](C)O")
#         result = [str(u) for u in sequence(mol)]
#         expected = ["B2", "A8", "A5", "AlphaAminoAcid", "B5", "B3", "C1", "B1", "D1"]
#         self.assertEqual(result, expected)

# class TestParsingTheopederinB(TestCase):
#     def test_linearization_theopederin_B(self):
#         mol = Molecule(r"C[C@H]1[C@H](O[C@](CC1=C)([C@H](C(=O)N[C@@H]2[C@@H]3[C@@H]([C@H](C([C@H](O3)C[C@@H](CCCC(=O)OC)O)(C)C)OC)OCO2)O)OC)C")
#         results = linearize(mol)
#         # Undirected ether linearization can lead to two different results.
#         expected1 = Molecule(r"C=C(CC(=O)[C@@H](O)C(=O)N[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(C)(C)C=C[C@H](O)CCCC(=O)O)[C@@H](C)[C@@H](C)O")
#         expected2 = Molecule(r"C=C(CC(=O)[C@@H](O)C(=O)NC(O)=C[C@H](O)[C@@H](O)C(C)(C)[C@H](O)C[C@H](O)CCCC(=O)O)[C@@H](C)[C@@H](C)O")
#         self.assertEqual(len(results), 2)
#         self.assertTrue(results[0] == expected1 or results[0] == expected2)
#         self.assertTrue(results[1] == expected1 or results[1] == expected2)

#     def test_sequencing_theopederin_B(self):
#         mol1 = Molecule(r"C=C(CC(=O)[C@@H](O)C(=O)N[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(C)(C)C=C[C@H](O)CCCC(=O)O)[C@@H](C)[C@@H](C)O")
#         mol2 = Molecule(r"C=C(CC(=O)[C@@H](O)C(=O)NC(O)=C[C@H](O)[C@@H](O)C(C)(C)[C@H](O)C[C@H](O)CCCC(=O)O)[C@@H](C)[C@@H](C)O")
#         result1 = [str(u) for u in sequence(mol1)]
#         result2 = [str(u) for u in sequence(mol2)]
#         expected1 = ["B2", "A8", "A5", "AlphaAminoAcid", "B5", "B3", "C1", "B1", "D1"]
#         expected2 = ["B2", "A8", "A5", "AlphaAminoAcid", "C", "B3", "B1", "B1", "D1"]
#         self.assertEqual(result1, expected1)
#         self.assertEqual(result2, expected2)

# class TestParsingThermolideA(TestCase):
#     def test_linearization_thermolide_A(self):
#         mol = Molecule(r"C[C@@H]1C[C@H]([C@@H](OC(=O)[C@H](NC(=O)C[C@H](C[C@@H]1O)O)C)[C@@H](C)C[C@H](C)[C@@H]([C@H](C)[C@@H](C[C@H](C)O)OC(=O)C)O)C")
#         results = linearize(mol)
#         expected = Molecule(r"C[C@@H]([C@@H](O)[C@@H](C)C[C@H](C)[C@H](O)[C@H](C)C[C@@H](C)[C@@H](O)C[C@H](O)CC(=O)N[C@H](C)C(=O)O)[C@H](O)C[C@H](C)O")
#         self.assertEqual(len(results), 1)
#         self.assertTrue(results[0] == expected)

#     def test_sequencing_thermolide_A(self):
#         mol = Molecule(r"C[C@@H]([C@@H](O)[C@@H](C)C[C@H](C)[C@H](O)[C@H](C)C[C@@H](C)[C@@H](O)C[C@H](O)CC(=O)N[C@H](C)C(=O)O)[C@H](O)C[C@H](C)O")
#         result = [str(u) for u in sequence(mol)]
#         expected = ["B1", "B2", "B2", "D2", "B2", "D2", "B1", "B1", "Ala"]
#         self.assertEqual(result, expected)

# class TestParsingThiocoraline(TestCase):
#     def test_linearization_thiocoraline(self):
#         mol = Molecule(r"CN1C2CSSCC(C(=O)N(C(C(=O)SCC(C(=O)NCC1=O)NC(=O)C3=NC4=CC=CC=C4C=C3O)CSC)C)N(C(=O)CNC(=O)C(CSC(=O)C(N(C2=O)C)CSC)NC(=O)C5=NC6=CC=CC=C6C=C5O)C")
#         results = linearize(mol)
#         expected = Molecule(r"O=C(CNC(=O)C(CS)NC(=O)c1nc2ccccc2cc1O)NC(CS)C(=O)NC(CS)C(=O)O")
#         self.assertEqual(len(results), 1)
#         self.assertTrue(results[0] == expected)

#     def test_sequencing_thiocoraline(self):
#         mol = Molecule(r"O=C(CNC(=O)C(CS)NC(=O)c1nc2ccccc2cc1O)NC(CS)C(=O)NC(CS)C(=O)O")
#         result = [str(u) for u in sequence(mol)]
#         expected = ["Cys", "Gly", "Cys", "Cys"]
#         self.assertEqual(result, expected)

# class TestParsingZincophorin(TestCase):
#     def test_linearization_zincophorin(self):
#         mol = Molecule(r"CCC[C@@H](C)/C=C(\C)/[C@@H]([C@H](C)/C=C/CC[C@H]([C@H](C)[C@@H]([C@H](C)[C@@H]([C@H](C)[C@@H]1[C@H](CC[C@H](O1)[C@H](C)C(=O)O)C)O)O)O)O")
#         results = linearize(mol)
#         # Undirected ether linearization can lead to two different results.
#         expected1 = Molecule(r"CCC[C@@H](C)/C=C(\C)[C@H](O)[C@H](C)/C=C/CC[C@@H](O)[C@H](C)[C@H](O)[C@H](C)[C@H](O)[C@H](C)[C@@H](O)[C@@H](C)CCC=C(C)C(=O)O")
#         expected2 = Molecule(r"CCC[C@@H](C)/C=C(\C)[C@H](O)[C@H](C)/C=C/CC[C@@H](O)[C@H](C)[C@H](O)[C@H](C)[C@H](O)C(C)=C[C@@H](C)CC[C@H](O)[C@H](C)C(=O)O")
#         self.assertEqual(len(results), 2)
#         self.assertTrue(results[0] == expected1 or results[0] == expected2)
#         self.assertTrue(results[1] == expected1 or results[1] == expected2)

#     def test_sequencing_zincophorin(self):
#         mol1 = Molecule(r"CCC[C@@H](C)/C=C(\C)[C@H](O)[C@H](C)/C=C/CC[C@@H](O)[C@H](C)[C@H](O)[C@H](C)[C@H](O)[C@H](C)[C@@H](O)[C@@H](C)CCC=C(C)C(=O)O")
#         mol2 = Molecule(r"CCC[C@@H](C)/C=C(\C)[C@H](O)[C@H](C)/C=C/CC[C@@H](O)[C@H](C)[C@H](O)[C@H](C)[C@H](O)C(C)=C[C@@H](C)CC[C@H](O)[C@H](C)C(=O)O")
#         result1 = [str(u) for u in sequence(mol1)]
#         result2 = [str(u) for u in sequence(mol2)]
#         expected1 = ["D2", "C2", "B2", "C1", "D1", "B2", "B2", "B2", "B2", "D1", "C2"]
#         expected2 = ["D2", "C2", "B2", "C1", "D1", "B2", "B2", "B2", "C2", "D1", "B2"]
#         self.assertEqual(result1, expected1)
#         self.assertEqual(result2, expected2)

# class TestParsingZwittermicinA(TestCase):
#     def test_linearization_zwittermicin_A(self):
#         mol = Molecule(r"C([C@H]([C@H](CO)N)O)[C@H]([C@H]([C@H]([C@@H](C(=O)N[C@@H](CNC(=O)N)C(=O)N)O)O)N)O")
#         results = linearize(mol)
#         expected = Molecule(r"NC[C@H](NC(=O)[C@@H](O)[C@H](O)[C@H](N)[C@H](O)C[C@@H](O)[C@@H](N)CO)C(N)=O")
#         self.assertEqual(len(results), 1)
#         self.assertTrue(results[0] == expected)

#     def test_sequencing_zwittermicin_A(self):
#         mol = Molecule(r"NC[C@H](NC(=O)[C@@H](O)[C@H](O)[C@H](N)[C@H](O)C[C@@H](O)[C@@H](N)CO)C(N)=O")
#         result = [str(u) for u in sequence(mol)]
#         expected = ["NRP", "B1", "B12", "B5", "AlphaAminoAcid"]
#         self.assertEqual(result, expected)