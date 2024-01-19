"""
Module for resolving primary sequence from RetroMol output.
"""
import typing as ty 

import networkx as nx 
import itertools 
from collections import defaultdict 
from copy import deepcopy 

from retromol.chem import ReactionTreeMapping, MonomerGraphMapping, MolecularPattern, identify_mol

def resolve_biosynthetic_sequence(
    reaction_tree: nx.DiGraph,
    reaction_mapping: ReactionTreeMapping, 
    monomer_graph: nx.Graph, 
    monomer_mapping: MonomerGraphMapping,
    motif_units: ty.List[MolecularPattern]
) -> ty.List[ty.Tuple[int, str]]:
    """
    Get depth-based biosynthetic sequence.
    
    :param nx.DiGraph reaction_tree: Reaction tree.
    :param ReactionTreeMapping reaction_mapping: Mapping of reaction products to molecules.
    :param nx.Graph monomer_graph: Monomer graph.
    :param MonomerGraphMapping monomer_mapping: Mapping of monomer graph nodes to molecules.
    :param ty.List[MolecularPattern] motif_units: List of molecular patterns.
    :returns: List of monomers in biosynthetic sequence.
    :rtype: ty.List[ty.Tuple[int, str]]
    """
    order = []

    for node in reaction_tree.nodes:
        if identity := identify_mol(reaction_mapping[node], motif_units):
            
            depth = 0
            current_node = node 
            while current_node in reaction_tree:

                current_node = list(reaction_tree.predecessors(current_node))
                if len(current_node) == 0: break 
                current_node = current_node[0]

                depth += 1
            
            order.append((node, depth, identity))

    order = sorted(order, key=lambda x: x[1], reverse=True)

    # Filter out nodes that are not in the monomer graph.
    order = [x for x in order if x[0] in monomer_mapping]

    # Convert reaction mappings to monomer mappings.
    order = [(monomer_mapping[x[0]][0], x[1], x[2]) for x in order]

    # Group by depth, as some nodes may have the same depth. This needs to be 
    # resolved.
    unresolved = defaultdict(list)
    for node, depth, identity in order:
        unresolved[depth].append((node, identity))
    
    resolved = []
    for depth in sorted(unresolved.keys(), reverse=True):

        # Seed resolved sequences. Every motif at lowest depth can be the first
        # motif in the sequence.
        if not len(resolved):
            resolved = [list(x) for x in itertools.permutations(unresolved[depth])]
            continue 

        # Find out which unit from the previous depth is the parent of the
        # current depth. 
        options = unresolved[depth]
        options = [x for x in itertools.permutations(options)]

        partially_resolved = []
        while resolved:

            seq = resolved.pop()
            
            # Check which unit is next by checking if extended sequence is a 
            # subgraph of the monomer graph.
            for option in options:
                option = [x for x in option]
                new_seq = deepcopy(seq)
                new_seq.extend(option)

                subgraph = nx.DiGraph()
                for i in range(len(new_seq) - 1):
                    subgraph.add_edge(new_seq[i][0], new_seq[i + 1][0])

                is_subgraph =  lambda G, H: all(h_edge in G.edges() for h_edge in H.edges)
                if is_subgraph(monomer_graph, subgraph):
                    partially_resolved.append(new_seq)

        if not len(partially_resolved):
            # No sequence could be extended. This is a dead end.
            return []
        
        resolved = partially_resolved

    if not len(resolved):
        return []
    else:
        return resolved[0]