import json
import logging
import typing as ty
from collections import defaultdict, deque
from copy import deepcopy

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import plotly.graph_objects as go
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

from retromol.retrosynthesis.parsing import parse_molecular_patterns
from retromol.retrosynthesis.chem import MolecularPattern, Molecule, ReactionRule
from retromol.retrosynthesis.graph import reaction_tree_to_monomer_graph, reaction_tree_to_digraph


monomers_path = r"app/src/server/data/monomers.json"
monomers_src = json.load(open(monomers_path, "r", encoding="utf-8"))
monomers = parse_molecular_patterns(monomers_src)


def mol_to_fingerprint(mol: Chem.Mol, radius: int, num_bits: int) -> np.array:
    fp_arr = np.zeros((0,), dtype=np.int8)
    fp_vec = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=num_bits)
    DataStructs.ConvertToNumpyArray(fp_vec, fp_arr)
    return fp_arr


def tanimoto_similarity(fp1: np.array, fp2: np.array) -> float:
    return np.logical_and(fp1, fp2).sum() / np.logical_or(fp1, fp2).sum()


def mol_to_encoding(mol: Chem.Mol, num_atoms: int, radius: int, num_bits: int) -> np.array:
    amns = [atom.GetIsotope() for atom in mol.GetAtoms() if atom.GetIsotope() > 0]
    amns = np.array([1 if x in amns else 0 for x in np.arange(num_atoms)])
    fp = mol_to_fingerprint(mol, radius, num_bits)
    return hash(np.hstack([fp, amns]).data.tobytes())


def find_roots(tree: defaultdict) -> ty.List[int]:
    roots = set(tree.keys())
    for children in tree.values():
        for products in children.values():
            for product in products:
                roots -= set(product) 
    return list(roots)


def find_leafs(tree: defaultdict) -> ty.List[int]:
    return [key for key, value in tree.items() if not value]


def count_edges(tree: defaultdict) -> int:
    return sum(len(children) for children in tree.values())


def visualize_tree(tree: defaultdict, mapping: dict) -> None:
    roots = find_roots(tree)
    leafs = find_leafs(tree)

    G = nx.DiGraph()

    for parent, children in tree.items():
        for reaction_idx, products in children.items():
            for product in products:
                for child in product:
                    G.add_edge(parent, child)

    pos = nx.spring_layout(G)
    for node in G.nodes(): G.nodes[node]['pos'] = list(pos[node])

    edge_x, edge_y = [], []
    for edge in G.edges():
        x0, y0 = G.nodes[edge[0]]['pos']
        x1, y1 = G.nodes[edge[1]]['pos']
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

    edge_trace = go.Scatter(x=edge_x, y=edge_y, line=dict(width=0.5, color='#888'), hoverinfo='none', mode='lines')

    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = G.nodes[node]['pos']
        node_x.append(x)
        node_y.append(y)

    color = ["red" if node in leafs else "green" if node in roots else "gray" for node in G.nodes]
    node_trace = go.Scatter(x=node_x, y=node_y, mode='markers', hoverinfo='text', marker=dict(showscale=False, color=color, size=10, line_width=2))
    
    node_text = []
    for node in G.nodes():
        mol = deepcopy(mapping[node])
        for atom in mol.GetAtoms(): atom.SetIsotope(0)
        smi = Chem.MolToSmiles(mol)
        node_text.append(smi)
    node_trace.text = node_text

    fig = go.Figure(data=[edge_trace, node_trace], layout=go.Layout(showlegend=False, hovermode='closest', margin=dict(b=20, l=5, r=5, t=40), xaxis=dict(showgrid=False, zeroline=False, showticklabels=False), yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))
    fig.show()


def apply_rules(mol: Chem.Mol, reactions: ty.List[str]) -> ty.Tuple[Chem.Mol, defaultdict, dict]:
    logger = logging.getLogger(__name__)

    rrs = [Chem.rdChemReactions.ReactionFromSmarts(rr) for rr in reactions]

    # sort rules
    # ... probably want to cycle between 1-to-1 and 1-to-many rules, start with 1-to-many apply until no more etc.

    radius = 2
    num_bits = 2048
    tree = defaultdict(lambda: defaultdict(set))
    mapping = dict()

    num_atoms = mol.GetNumAtoms()
    for atom in mol.GetAtoms():
        atom.SetIsotope(atom.GetIdx() + 1)

    reaction_matches = set()

    mols = [mol]
    while mols:
        current = mols.pop()
        current_encoding = mol_to_encoding(current, num_atoms, radius, num_bits)
        tree[current_encoding]
        mapping[current_encoding] = current

        for rr_idx, rr in enumerate(rrs):
            
            # Get the reactant part of the reaction.
            matching_part = Chem.MolFromSmarts(reactions[rr_idx].split(">>")[0])
            matches = current.GetSubstructMatches(matching_part)
            
            # Check if the current molecule matches the reactant reaction pattern.
            if len(matches) > 0:
                reaction_matches.add(rr_idx)
            else:
                continue

            # If reaction matches are ambiguous, we greedily select the first ones.
            reaction_results = rr.RunReactants([current])
            if len(matches) != len(reaction_results):
                reaction_results = reaction_results[:len(matches)]

            for results in reaction_results:

                if len(results) > 0:
                    reaction_products = list()

                    # All the results are sanitized.
                    sanitized_results = list()
                    for result in results:
                        try:
                            Chem.SanitizeMol(result)
                            sanitized_results.append(result)
                        except Exception:
                            msg = f"Failed to sanitize the molecule: {Chem.MolToSmiles(result)}"
                            logger.error(msg)
                            continue  # Quick fix for sanitization issues.
                    
                    # Only add the results that were fully successfully sanitized.
                    if len(results) == len(sanitized_results):
                        for result in sanitized_results:
                            result_encoding = mol_to_encoding(result, num_atoms, radius, num_bits)
                            reaction_products.append(result_encoding)

                            if result_encoding not in tree:
                                mols.append(result)

                    tree[current_encoding][rr_idx].add(frozenset(reaction_products))

    logger.debug(f"Tried to apply {len(reaction_matches)} reaction rules to the molecule:")
    for reaction_name in reaction_matches:
        logger.debug(f" >> {reaction_name}")

    return tree, mapping


# def invert_tree(tree: defaultdict) -> defaultdict:
#     inverted = defaultdict(lambda: defaultdict(set))
#     for parent, children in tree.items():
#         for reaction, products in children.items():
#             for product in products:
#                 for child in product:
#                     inverted[child][reaction].add(parent)

#     return inverted

def get_routes(tree: defaultdict, node_a: int, node_b: int) -> ty.List[ty.List[int]]:
    # tree = reaction_tree_to_digraph(tree)

    def find_all_paths(graph, start_node, target_node):
        def dfs(current_node, target_node, path, all_paths, visited):
            path.append(current_node)
            visited.add(current_node)

            if current_node == target_node:
                all_paths.append(list(path))
            else:
                # for neighbor in graph.successors(current_node):
                for reaction in graph[current_node]:
                    neighbors = graph[current_node][reaction]
                    for neighbor_set in neighbors:
                        for neighbor in neighbor_set:
                            if neighbor not in visited:  # Avoid cycles
                                dfs(neighbor, target_node, path, all_paths, visited)

            # Backtrack
            path.pop()
            visited.remove(current_node)

        all_paths = []
        dfs(start_node, target_node, [], all_paths, set())
        return all_paths
    
    return find_all_paths(tree, node_a, node_b)


def find_shortest_paths(reaction_tree, root, targets, mapping):
    """
    Find the shortest paths to the target nodes in the reaction tree.

    Parameters:
    - reaction_tree: dict, where keys are reactants and values are dicts mapping reaction indices to sets of frozensets of products.
    - root: The root reactant from which to start the search.
    - targets: A set of target reactants you want to reach.

    Returns:
    A dict mapping each target to the shortest sequence of reaction indices that produces it.
    """
    # Initialize a queue for BFS
    queue = deque([(root, [])])
    # Initialize a dictionary to store the shortest paths
    shortest_paths = {}

    while queue:
        current_node, path = queue.popleft()

        # If the current node is a target, store the path and stop searching for this target
        if current_node in targets:
            shortest_paths[current_node] = path
            # Remove the target from the set to stop looking for it
            targets.remove(current_node)
            # If we've found all targets, we can return the results
            if not targets:
                return shortest_paths

        # Get the reactions that can be applied to the current node
        for reaction_index, products in reaction_tree.get(current_node, {}).items():
            # For each product set, consider each product (there could be multiple in a frozenset)
            for product_set in products:

                for product in product_set:
                    # Add the product to the BFS queue with the updated path
                    queue.append((product, path + [reaction_index]))

    # Return the shortest paths found
    return shortest_paths


def find_shortest_paths2(reaction_tree, root, targets, mapping):
    """
    Find the shortest path that reaches all target nodes in the reaction tree.

    Parameters:
    - reaction_tree: dict, where keys are reactants and values are dicts mapping reaction indices to sets of frozensets of products.
    - root: The root reactant from which to start the search.
    - targets: A set of target reactants you want to reach.

    Returns:
    A list representing the shortest sequence of reaction indices that produces all targets.
    """
    # Initialize a queue for BFS
    queue = deque([(root, [])])
    # Initialize a set to store the visited nodes
    visited = set()

    while queue:
        current_node, path = queue.popleft()

        # If the current node is a target, add it to the visited set
        if current_node in targets:
            visited.add(current_node)
            # If we've found all targets, we can return the path
            if visited == targets:
                return path

        # Get the reactions that can be applied to the current node
        for reaction_index, products in reaction_tree.get(current_node, {}).items():
            # For each product set, consider each product (there could be multiple in a frozenset)
            for product_set in products:
                for product in product_set:
                    # Add the product to the BFS queue with the updated path
                    queue.append((product, path + [reaction_index]))

    # If no path includes all targets, return an empty list
    return []


def prune_tree(tree: defaultdict, mapping: dict) -> None:
    roots = find_roots(tree)
    if len(roots) != 1: raise ValueError("Tree must have exactly one root.")
    root = roots[0]
    mol = Molecule("test", "*")
    mol.smiles = Chem.MolToSmiles(mapping[root])
    mol.compiled = mapping[roots[0]]
    monomer_graph, monomer_graph_mapping = reaction_tree_to_monomer_graph(mol, reaction_tree_to_digraph(tree), mapping, monomers, return_unmatched=True) 
    mapped_nodes = monomer_graph_mapping.keys()

    # Get the indices of the mapped atoms in each molecule
    mapping_inds = {}
    for k in mapping:
        mol = deepcopy(mapping[k])
        inds = []
        for atom in mol.GetAtoms():
            if atom.GetIsotope() > 0:
                inds.append(atom.GetIsotope())
        mapping_inds[k] = inds

    targets = set(mapped_nodes)

    # remove redundant nodes
    to_remove = set()
    for node in tree.keys():
        if node not in targets:
            x = deepcopy(mapping[node])
            inds = mapping_inds[node]
            # if inds is a subset of a target, remove the node
            for target in targets:
                if set(inds).issubset(mapping_inds[target]):
                    to_remove.add(node)
                    break

    new_tree = defaultdict(lambda: defaultdict(set))

    for node in tree.keys():
        if node not in to_remove:
            for reaction, products in tree[node].items():
                for product in products:
                    if all(p not in to_remove for p in product):
                        new_tree[node][reaction].add(product)

    # Find the shortest paths to the target nodes
    shortest_paths = find_shortest_paths(new_tree, root, targets, mapping)
    
    for k, reaction_indices in shortest_paths.items():
        x = deepcopy(mapping[k])
        for atom in x.GetAtoms(): atom.SetIsotope(0)
        print(k, reaction_indices, Chem.MolToSmiles(x))

    return new_tree
    


def main() -> None:
    smi1 = r"CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O"
    smi2 = r"CCCCCCCCCC(=O)NC(CC1=CNC2=CC=CC=C21)C(=O)NC(CC(=O)N)C(=O)NC(CC(=O)O)C(=O)NC3C(OC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)C(C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C"

    # for smi2 (daptomycin): shouldn't be able to use reaction routes that result in contradictions (products that have overlapping atoms but are not identified from greedy set cover)

    rr1 = r"[C,c:1][C;!R:2]~[C;!R:3][C:4](=[O:5])[OH:6]>>[C:1]C(=O)[OH].[OH][S][C:2]~[C:3][C:4](=[O:5])[OH:6]" # pks
    rr2 = r"[C,c;R:1]-[C;R:2](=[O:3])-[O;R:4]-[C,c;R:5]>>([C,c:1]-[C:2](=[O:3])-[OH].[OH:4]-[C,c:5])" # macrocyclization 
    rr3 = r"[C:1][O:2][C:3]1[O:4][C:5][C:6][C:7][C:8]1>>[C:1][OH:2].[OH][C:3]1[O:4][C:5][C:6][C:7][C:8]1" # deglycosylation
    rr4 = r"[C:1][O:2][CH3:3]>>[C:1][OH:2].[CH4:3]" # demethylation
    rr5 = r"[C:1][NH1:2][CH3:3]>>[C:1][NH2:2].[CH4:3]" # demethylation
    rr6 = r"[C:1][NH0:2]([C:4])[CH3:3]>>[C:1][NH1:2][C:4].[CH4:3]" # demethylation
    rr7 = r"[*:1][C:2](=[O:3])[NH1:4][C:5][C:6](=[O:7])[OH:8]>>[C:1][C:2](=[O:3])[OH].[NH2:4][C:5][C:6](=[O:7])[OH:8]"
    rr8 = r"[*:1][C:2](=[O:3])[NH0:4][C:5][C:6](=[O:7])[OH:8]>>[C:1][C:2](=[O:3])[OH].[NH1:4][C:5][C:6](=[O:7])[OH:8]"
    reactions = [rr1, rr2, rr3, rr4, rr5, rr6, rr7, rr8]

    print("\nerythromycin:")
    tree, mapping = apply_rules(Chem.MolFromSmiles(smi1), reactions)
    pruned = prune_tree(tree, mapping)

    print("\ndaptomycin:")
    tree, mapping = apply_rules(Chem.MolFromSmiles(smi2), reactions)
    pruned = prune_tree(tree, mapping)

    exit(0)


if __name__ == "__main__":
    main()
