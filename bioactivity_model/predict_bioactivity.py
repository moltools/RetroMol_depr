import joblib 

import numpy as np
from rdkit import Chem 

from retromol.chem import mol_to_fingerprint, draw_molecule

radius = 4
num_bits = 2048

smi_complete = r"CCC1OC(=O)C(C)C(OC2CC(C)(OC)C(O)C(C)O2)C(C)C(OC2OC(C)CC(N(C)C)C2O)C(C)(O)CC(C)C(=O)C(C)C(O)C1(C)O"
mol_complete = Chem.MolFromSmiles(smi_complete)
fp_complete = mol_to_fingerprint(mol_complete, radius, num_bits)

smi_no_sugars = r"CCC1OC(=O)C(C)C(O)C(C)C(O)C(C)(O)CC(C)=C(O)C(C)C(O)C1(C)O"
mol_no_sugars = Chem.MolFromSmiles(smi_no_sugars)
fp_no_sugars = mol_to_fingerprint(mol_no_sugars, radius, num_bits)

smi_design = r"CCC1OC(=O)C(C)C(OC2OC(CO)C(O)C(O)C(O)2)C(C)C(O)C(C)(O)CC(C)=C(O)C(C)C(O)C1(C)O" # beter
# smi_design = r"CCC1OC(=O)C(C)C(O)C(C)C(O)C(C)(O)CC(C)=C(OC2OC(CO)C(O)C(O)C(O)2)C(C)C(O)C1(C)O" # slechter
mol_design = Chem.MolFromSmiles(smi_design)
fp_design = mol_to_fingerprint(mol_design, radius, num_bits)

model = joblib.load('model.joblib')

pred_complete = model.predict_proba(np.array([fp_complete]))[0][1]
print(pred_complete)

pred_no_sugars = model.predict_proba(np.array([fp_no_sugars]))[0][1]
print(pred_no_sugars)

pred_design = model.predict_proba(np.array([fp_design]))[0][1]
print(pred_design)

draw_molecule(mol_complete, f"erythromycin_complete_{pred_complete}.png")
draw_molecule(mol_no_sugars, f"erythromycin_no_sugars_{pred_no_sugars}.png")
draw_molecule(mol_design, f"erythromycin_design_{pred_design}.png")