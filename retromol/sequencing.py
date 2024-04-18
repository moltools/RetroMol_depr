"""This module contains functions for retrieving monomer sequences from 
monomer graphs.
"""
import typing as ty

import networkx as nx

from retromol.parsing import Result

def parse_modular_natural_product(rec: Result) -> ty.List[ty.List[ty.Any]]:
    """Parse a monomer graph to retrieve a modular natural product sequence.
    
    :param rec: The RetroMol Result object.
    :type rec: Result
    :return: The modular natural product sequence.
    :rtype: ty.List[ty.List[ty.Any]]
    """
    # Map reaction_tree_id to monomer_id.
    monomer_mapping = {}
    for monomer_id, items in rec.monomer_graph.items():
        if items["identity"] is not None:
            monomer_mapping[items["reaction_tree_id"]] = monomer_id

    # Turn the monomer graph into an undirected graph.
    monomer_graph = {}
    for _, items in rec.monomer_graph.items():
        if items["identity"] is not None:
            node = items["reaction_tree_id"]
            monomer_graph[node] = []

            for neighbor in items["neighbors"]:
                if rec.monomer_graph[neighbor]["identity"] is not None:
                    neighbor_node = rec.monomer_graph[neighbor]["reaction_tree_id"]
                    monomer_graph[node].append(neighbor_node)

    # Turn reaction tree into digraph.
    reaction_tree = nx.DiGraph()
    for node, items in rec.reaction_tree.items():
        children = items["children"]
        reaction_tree.add_node(node)
        for child in children:
            reaction_tree.add_edge(node, child)

    # Get root from reaction tree.
    root = [node for node, degree in reaction_tree.in_degree() if degree == 0][0]

    # Get depths for all identified monomers.
    monomer_depths = {}
    for node in monomer_graph:
        depth = nx.shortest_path_length(reaction_tree, source=root, target=node)
        monomer_depths[node] = depth

    # Get shallowest node.
    shallowest = min(monomer_depths, key=monomer_depths.get)

    # Find non-overlapping path from shallowest to all other nodes.
    monomer_graph = nx.Graph(monomer_graph)
    path = [x for x in nx.dfs_preorder_nodes(monomer_graph, shallowest)]
    path = path[::-1]

    # Check if there are left-over nodes.
    if len(path) != len(monomer_graph):
        return [] # TODO: implement this.

    else:
        monomer_ids = [monomer_mapping[node] for node in path]
        monomer_props = [rec.monomer_graph[monomer_id]["identity"] for monomer_id in monomer_ids]
        return [monomer_props]
