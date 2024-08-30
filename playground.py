import json
import itertools
import logging
import typing as ty
from collections import defaultdict, deque
from copy import deepcopy

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import plotly.graph_objects as go
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize

from retromol.retrosynthesis.parsing import parse_molecular_patterns, parse_reaction_rules
from retromol.retrosynthesis.chem import MolecularPattern, Molecule, ReactionRule
from retromol.retrosynthesis.graph import reaction_tree_to_monomer_graph, reaction_tree_to_digraph


RDLogger.DisableLog('rdApp.*')


monomers_path = r"app/src/server/data/monomers.json"
monomers_src = json.load(open(monomers_path, "r", encoding="utf-8"))
monomers = parse_molecular_patterns(monomers_src)

reactions_path = r"app/src/server/data/reactions.json"
reactions_src = json.load(open(reactions_path, "r", encoding="utf-8"))
reactions = parse_reaction_rules(reactions_src)


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

    return mol, tree, mapping

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


def find_shortest_paths(reaction_tree, root, targets, mapping, trans_map):
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
                    
                    # unique identifier for input->reaction->output
                    key = (current_node, reaction_index, product_set)
                    reaction_identifier = trans_map[key]

                    # Add the product to the BFS queue with the updated path
                    queue.append((product, path + [(reaction_index, reaction_identifier)]))

    # Return the shortest paths found
    return shortest_paths


def prune_tree(tree: defaultdict, mapping: dict, trans_map: dict) -> None:
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
    shortest_paths = find_shortest_paths(new_tree, root, deepcopy(targets), mapping, trans_map)
    shortest_paths = [(k, v[::-1]) for k, v in shortest_paths.items()]

    # for k, v in shortest_paths:
    #     print(k, v)

    # print(targets)
    # for target in targets:
    #     print(target)
    #     all_paths = get_routes(new_tree, root, target)
    #     print(all_paths)

    # exit()

    # sort from longest to shortest path
    shortest_paths = sorted(shortest_paths, key=lambda x: len(x[1]), reverse=True)

    # Identify the monomers
    identified = {}
    for k, v in monomer_graph_mapping.items():
        identified[k] = monomer_graph_mapping[k][1]

    return new_tree, shortest_paths, identified
    

def reverse_reaction(rxn_repr: str):
    rxn_repr = rxn_repr.split(">>")
    new_repr = rxn_repr[1] + ">>" + rxn_repr[0]
    compiled = Chem.rdChemReactions.ReactionFromSmarts(new_repr)
    return compiled


def add_path(graph, path):
    node = graph
    for value in path:
        if value not in node:
            node[value] = {}
        node = node[value]

def merge_paths(paths):
    if not paths:
        return {}

    graph = {}
    for path in paths:
        node = graph
        for value in path:
            if value in node:
                node = node[value]
            else:
                break
        else:
            continue
        add_path(node, path[len(path) - len(node):])

    return graph

def is_subsequence(a, b):
    if len(a) == 0 or len(b) == 0:
        return False
    
    if len(a) == len(b):
        return False 

    iter_a = iter(a)
    return all(item in iter_a for item in b)

def reconstruct_synthesis(mol, synth, pruned, mapping, rev_reactions, identified):    
    building_blocks = []
    for x, _ in synth:
        y = deepcopy(mapping[x])
        # for atom in mol.GetAtoms():
        #     atom.SetIsotope(0)
        building_blocks.append((x, y))

    fps = {}
    for n in pruned.keys():
        x = deepcopy(mapping[n])
        for atom in x.GetAtoms():
            atom.SetIsotope(0)
        fps[n] = mol_to_fingerprint(x, 2, 2048)

    graph_nodes = []

    print(building_blocks)
    for _, path in synth:
        for rxn_idx in path:
            rxn = rev_reactions[rxn_idx[0]]
            num_inputs = len(rxn.GetReactants())
            perms = itertools.permutations(building_blocks, num_inputs)
            perms = [p for p in perms if len(p) == num_inputs]
            for perm in perms:
                out = rxn.RunReactants([p[1] for p in perm])
                for result in out:
                    # result needs to be 1 output
                    if len(result) != 1: raise ValueError("Reaction must have exactly one output.")
                    result = result[0]

                    # sanitize, but ignore valence
                    for atom in result.GetAtoms():
                        atom.SetNoImplicit(False)
                        atom.SetNumExplicitHs(0)
                        # atom.SetIsotope(0)
                    # Chem.SanitizeMol(result)
                    Chem.AssignStereochemistry(result, cleanIt=True)

                    try:
                        Chem.SanitizeMol(result)
                    except Exception:
                        continue
                    # fp = mol_to_fingerprint(result, 2, 2048)

                    enc = mol_to_encoding(result, mol.GetNumAtoms(), 2, 2048)
                    if enc in pruned:
                        # removed used building blocks
                        used = [p[0] for p in perm]
                        building_blocks = [p for p in building_blocks if p[0] not in used]

                        # add result as new building block
                        building_blocks.append((enc, result))

                        # print smiles of result
                        q = deepcopy(result)
                        for atom in q.GetAtoms():
                            atom.SetIsotope(0)
                        graph_nodes.append(q)
                        smi = Chem.MolToSmiles(q)
                        print(smi)

                        break
        
        break


    print(graph_nodes)  

    # calculate reconstruction error (tanimoto distance final graph node and input)
    fp_input = mol_to_fingerprint(mol, 2, 2048)
    fp_reconstructed = mol_to_fingerprint(graph_nodes[-1], 2, 2048)
    print("\n\nRECONSTRUCTION ERROR:", tanimoto_similarity(fp_input, fp_reconstructed))
    print("\n\n")

    # every graph node is molecule, draw as image and add to graph from left to right
    # add edges between nodes (directed)

    # g = nx.DiGraph()
    # for i in range(len(graph_nodes)):
    #     g.add_node(i, smiles=Chem.MolToSmiles(graph_nodes[i]))

    # for i in range(len(graph_nodes) - 1):
    #     g.add_edge(i, i + 1)

    # pos = nx.spring_layout(g)
    # for node in g.nodes(): g.nodes[node]['pos'] = list(pos[node])
    
    # edge_x, edge_y = [], []
    # for edge in g.edges():
    #     x0, y0 = g.nodes[edge[0]]['pos']
    #     x1, y1 = g.nodes[edge[1]]['pos']
    #     edge_x.append(x0)
    #     edge_x.append(x1)
    #     edge_x.append(None)
    #     edge_y.append(y0)
    #     edge_y.append(y1)
    #     edge_y.append(None)

    # edge_trace = go.Scatter(x=edge_x, y=edge_y, line=dict(width=0.5, color='#888'), hoverinfo='none', mode='lines')

    # node_x = []
    # node_y = []
    # for node in g.nodes():
    #     x, y = g.nodes[node]['pos']
    #     node_x.append(x)
    #     node_y.append(y)

    # color = ["red" if node in identified else "green" for node in g.nodes]
    # node_trace = go.Scatter(x=node_x, y=node_y, mode='markers', hoverinfo='text', marker=dict(showscale=False, color=color, size=10, line_width=2))


    # node_text = []
    # for node in g.nodes():
    #     node_text.append(g.nodes[node]['smiles'])
    # node_trace.text = node_text

    # fig = go.Figure(data=[edge_trace, node_trace], layout=go.Layout(showlegend=False, hovermode='closest', margin=dict(b=20, l=5, r=5, t=40), xaxis=dict(showgrid=False, zeroline=False, showticklabels=False), yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))
    # fig.show()


def merge_paths(synthesis):
    paths = [p[1] for p in synthesis]
    merged = "->".join([f"{x[0]}+{x[1]}" for x in paths[0]])
    others = []
    for path in paths[1:]:
        others.append("->".join([f"{x[0]}+{x[1]}" for x in path]))
    
    unmerged = []
    for i in range(len(others)):
        other = others[i]
        if other in merged:
            continue
        else:
            unmerged.append(other)


    print(merged)
    print(unmerged)
    print("\n\n")


def main() -> None:
    smi1 = r"CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O"
    smi2 = r"CCCCCCCCCC(=O)NC(CC1=CNC2=CC=CC=C21)C(=O)NC(CC(=O)N)C(=O)NC(CC(=O)O)C(=O)NC3C(OC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)C(C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C"
    smi3 = r"CC1C(OC(C(C1OC)(C)C)(C(C)CC(C)C(C2(C(O2)C(C)C=C(C)C)C)O)O)CC(=O)O"

    rr1 = r"[C,c:1][C;!R:2]~[C;!R:3][C:4](=[O:5])[OH:6]>>[C:1]C(=O)[OH].[OH][S][C:2]~[C:3][C:4](=[O:5])[OH:6]" # pks
    rr2 = r"[C,c;R:1]-[C;R:2](=[O:3])-[O;R:4]-[C,c;R:5]>>([C,c:1]-[C:2](=[O:3])-[OH].[OH:4]-[C,c:5])" # macrocyclization 
    rr3 = r"[C:1][O:2][C:3]1[O:4][C:5][C:6][C:7][C:8]1>>[C:1][OH:2].[OH][C:3]1[O:4][C:5][C:6][C:7][C:8]1" # deglycosylation
    rr4 = r"[C:1][O:2][CH3:3]>>[C:1][OH:2].[CH4:3]" # demethylation
    rr5 = r"[C:1][NH1:2][CH3:3]>>[C:1][NH2:2].[CH4:3]" # demethylation
    rr6 = r"[C:1][NH0:2]([C:4])[CH3:3]>>[C:1][NH1:2][C:4].[CH4:3]" # demethylation
    rr7 = r"[*:1][C:2](=[O:3])[NH1:4][C:5][C:6](=[O:7])[OH:8]>>[C:1][C:2](=[O:3])[OH].[NH2:4][C:5][C:6](=[O:7])[OH:8]"
    rr8 = r"[*:1][C:2](=[O:3])[NH0:4][C:5][C:6](=[O:7])[OH:8]>>[C:1][C:2](=[O:3])[OH].[NH1:4][C:5][C:6](=[O:7])[OH:8]"
    rr9 = r"[*:1]-[CH0:2]1(-[OH:3])-[O:4]-[CH1:5](-[*:6])-[C:7][C:8][C:9]1>>([*:1]-[CH0:2]1=[OH0:3].[OH1:4]-[CH1:5](-[*:6])-[C:7][C:8][C:9]1)" # ether
    rr10 = r"[C:1]1[O:2][C:3]1>>[C:1]=[C:3].[O:2]" # epoxidation 
    reactions = [rr1, rr2, rr3, rr4, rr5, rr6, rr7, rr8, rr9, rr10]
    
    mol, tree, mapping = apply_rules(Chem.MolFromSmiles(smi2), reactions)


    # map every inputKey-reaciton-outputKey to a unique identifier
    transformation_mapping = {}
    for inputKey in tree:
        for reaction, outputKeys in tree[inputKey].items():
            for outputKey in outputKeys:
                transformation_mapping[(inputKey, reaction, outputKey)] = len(transformation_mapping)

    pruned, synthesis, identified = prune_tree(tree, mapping, transformation_mapping)

    print("\n")
    for n, s in synthesis:
        print(n, s)
    print("\n")

    # NOTE:
    # this merging with reaction nodes is good approach, just make sure you get more
    # than just shortest paths. Need all paths and then greedy approach for reconstruction
    # ideally you would see with multi-mers or compounds with multiple parts that the merged
    # pipeline would split somewhere
    merge_paths(synthesis)

    exit()

    rev_reactions = [reverse_reaction(rr) for rr in reactions]
    forward_synthesis = reconstruct_synthesis(mol, synthesis, pruned, mapping, rev_reactions, identified)
    print(forward_synthesis)

    exit(0)


if __name__ == "__main__":
    main()
