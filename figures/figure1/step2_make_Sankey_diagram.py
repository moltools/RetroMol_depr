"""
Make a Sankey diagram.

Dependencies:
    - matplotlib
    - sankeyflow
"""
import matplotlib.pyplot as plt
from sankeyflow import Sankey

num_compounds_npatlas = 33372
num_compounds_selected = 8276
num_compounds_not_selected = num_compounds_npatlas - num_compounds_selected

num_compounds_successfully_parsed = 8023 
num_compounds_failed_to_parse = 254

num_compounds_fully_identified = 766
num_compounds_fully_unidentified = 3068
num_compounds_partially_identified = num_compounds_successfully_parsed - num_compounds_fully_identified - num_compounds_fully_unidentified

# Make Sankey diagram.
nodes = [
    [("NPAtlas", num_compounds_npatlas)],
    [("Modular NPs", num_compounds_selected), ("Other", num_compounds_not_selected)],
    [("Parsed", num_compounds_successfully_parsed), ("Failed", num_compounds_failed_to_parse)],
    [("Fully identified", num_compounds_fully_identified), ("Partially identified", num_compounds_partially_identified), ("Fully unidentified", num_compounds_fully_unidentified)],
]
flows = [
    ("NPAtlas", "Modular NPs", num_compounds_selected),
    ("NPAtlas", "Other", num_compounds_not_selected),
    ("Modular NPs", "Parsed", num_compounds_successfully_parsed),  
    ("Modular NPs", "Failed", num_compounds_failed_to_parse),
    ("Parsed", "Fully identified", num_compounds_fully_identified),
    ("Parsed", "Partially identified", num_compounds_partially_identified),
    ("Parsed", "Fully unidentified", num_compounds_fully_unidentified),
] 

plt.figure(figsize=(4, 3), dpi=144)
s = Sankey(flows=flows, nodes=nodes)
s.draw()
plt.tight_layout()
plt.show()