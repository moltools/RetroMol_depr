import itertools
import typing as ty 
from collections import defaultdict 
from copy import deepcopy 

import networkx as nx 

from retromol.chem import (
    Tree, ReactionTreeMapping, MonomerGraphMapping, Molecule, MolecularPattern, 
    identify_mol,
    greedy_max_set_cover
)

def reaction_tree_to_digraph(tree: Tree) -> nx.DiGraph:
    """
    Convert reaction tree to directed graph.

    Parameters
    ----------
    tree : Tree
        Reaction tree.
    
    Returns
    -------
    digraph : nx.DiGraph
        Directed graph.
    """
    digraph = nx.DiGraph()

    for parent, reactions in tree.items():
        digraph.add_node(parent)

        for reaction, children in reactions.items():
            for child in children:
                for product in child:
                    digraph.add_node(product)
                    digraph.add_edge(parent, product, reaction=reaction)

    return digraph

def reaction_tree_to_monomer_graph(
    mol: Molecule, 
    tree: nx.DiGraph,
    mapping: ReactionTreeMapping,
    monomers: ty.List[MolecularPattern]
) -> ty.Tuple[nx.Graph, MonomerGraphMapping]:
    """
    Convert reaction tree to monomer graph.

    Parameters
    ----------
    mol : Molecule
        Molecule.
    tree : nx.DiGraph
        Reaction tree.
    mapping : ReactionTreeMapping
        Mapping.
    monomers : ty.List[MolecularPattern]
        List of monomers.
    
    Returns
    -------
    monomer_graph : nx.Graph
        Monomer graph.
    monomer_graph_mapping : MonomerGraphMapping
        Monomer graph mapping.
    """
    identified = []
    for node in tree.nodes:
        if identity := identify_mol(mapping[node], monomers):
            identified.append((node, identity))
    
    subgraphs = greedy_max_set_cover(mol.compiled, identified, mapping)

    # All atoms in the substrate have a valid atom map number as isotope.
    monomer_graph = nx.Graph()
    for atom in mol.compiled.GetAtoms():
        monomer_graph.add_node(atom.GetIsotope())
    for bond in mol.compiled.GetBonds():
        monomer_graph.add_edge(
            bond.GetBeginAtom().GetIsotope(),
            bond.GetEndAtom().GetIsotope()
        )
    
    monomer_graph_mapping = dict()
    for subgraph in subgraphs:
        amns = [
            atom.GetIsotope()
            for atom in mapping[subgraph[0]].GetAtoms()
            if atom.GetIsotope() > 0
        ]

        # Merge subgraph nodes in monomer graph.
        monomer_graph.add_node(amns[0])
        for amn in amns[1:]:
            monomer_graph = nx.contracted_nodes(
                monomer_graph,
                amns[0],
                amn,
                self_loops=False
            )
        
        monomer_graph_mapping[subgraph[0]] = (amns[0], subgraph[1])

    return monomer_graph, monomer_graph_mapping 

def resolve_biosynthetic_sequence(
    reaction_tree: nx.DiGraph,
    reaction_mapping: ReactionTreeMapping, 
    monomer_graph: nx.Graph, 
    monomer_mapping: MonomerGraphMapping,
    motif_units: ty.List[MolecularPattern]
) -> ty.List[ty.Tuple[int, str]]:
    """
    Get depth-based biosynthetic sequence.
    
    Parameters
    ----------
    reaction_tree : nx.DiGraph
        Reaction tree.
    reaction_mapping : ReactionTreeMapping
        Reaction mapping.
    monomer_graph : nx.Graph
        Monomer graph.
    monomer_mapping : MonomerGraphMapping
        Monomer mapping.
    motif_units : ty.List[MolecularPattern]
        List of motif units.
    
    Returns
    -------
    order : ty.List[ty.Tuple[int, str]]
        List of tuples containing node id in monomer graph and identity.
    """
    # # Monomer graph node id to reaction node id.
    # monomer_to_reaction = {v[0]: k for k, v in monomer_mapping.items()}

    # # Check if node is a core node.
    # def is_core_node(name: str) -> bool:
    #     for unit in motif_units:
    #         if name == unit.name: return True
    #     return False

    # # Core nodes. 
    # core_nodes = [
    #     monomer_id for _, (monomer_id, identity) 
    #     in monomer_mapping.items() 
    #     if is_core_node(identity)
    # ]

    # # Create a monomer graph with only core nodes.
    # core_monomer_graph = nx.Graph()
    # for node in core_nodes:
    #     core_monomer_graph.add_node(node)
    # for edge in monomer_graph.edges:
    #     if edge[0] in core_nodes and edge[1] in core_nodes:
    #         core_monomer_graph.add_edge(edge[0], edge[1])
    
    # # Identify all possible paths from core nodes to core nodes.
    # paths = []
    # for node in core_nodes:
    #     for other_node in core_nodes:
    #         if node == other_node: continue 
    #         for path in nx.all_simple_paths(core_monomer_graph, node, other_node):
    #             if len(path) == len(core_monomer_graph.nodes):
    #                 paths.append(path)

    # # Calculate depth of each monomer graph node in the reaction tree.
    # depths = dict()
    # for node in core_monomer_graph.nodes:
    #     node_id = monomer_to_reaction[node]
    #     depth = 0
    #     current_node = node_id
    #     while current_node in reaction_tree:
    #         current_node = list(reaction_tree.predecessors(current_node))
    #         if len(current_node) == 0: break 
    #         current_node = current_node[0]
    #         depth += 1
    #     depths[node] = depth

    # # Get the path that has the longest sequence of monomers where depth 
    # # increases every step.
    # path = None
    # for p in paths:
    #     if path is None:
    #         path = p 
    #         continue 
    #     if len(p) > len(path):
    #         path = p 
    #         continue 
    #     if len(p) == len(path):
    #         if all(depths[p[i]] > depths[path[i]] for i in range(len(p))):
    #             path = p

    # # Get identity of each monomer graph node.
    # identities = dict()
    # for node in core_monomer_graph.nodes:
    #     identities[node] = monomer_mapping[monomer_to_reaction[node]][1]
    
    # # Get identities of each monomer graph node in the path.
    # path = [identities[node] for node in path]

    # print(path)

    # return path

    # # # Get all possible sequences.
    # # print(path)
    
    # # return []
    

    # # ====

    # Get the biosynthetic sequence solely based on depth of motifs in the
    # reaction tree.
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

