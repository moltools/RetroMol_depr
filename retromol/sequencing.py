"""This module contains functions for retrieving monomer sequences from
monomer graphs.
"""

import logging
import typing as ty

import networkx as nx
from rdkit import Chem


def parse_modular_natural_product(
    reaction_tree: ty.Dict[int, ty.List[int]],
    monomer_graph: ty.Dict[str, ty.Any],
    core_types: ty.List[str] = ["polyketide", "peptide"]
) -> ty.List[ty.Tuple[Chem.Mol, ty.List[str]]]:
    """Parse a monomer graph to retrieve a modular natural product sequence.

    :param reaction_tree: The reaction tree.
    :type reaction_tree: ty.Dict[int, ty.List[int]]
    :param monomer_graph: The monomer graph.
    :type monomer_graph: ty.Dict[str, ty.Any]
    :param core_types: The core types to consider.
    :type core_types: ty.List[str]
    :return: The modular natural product sequences.
    :rtype: ty.List[ty.Tuple[Chem.Mol, ty.List[str]]]
    """
    logger = logging.getLogger(__name__)

    # Get all identified monomers as Chem.Mol objects.
    monomer_ids = []
    for _, props in monomer_graph.items():
        identity = props["identity"]
        if identity is not None:
            if identity.split("|")[0] not in core_types:
                continue
            monomer_ids.append((props["reaction_tree_id"], identity))
    subs = [(Chem.MolFromSmiles(reaction_tree[x]["smiles"]), identity) for x, identity in monomer_ids]

    # Retrieve the original AMNs of the identified core type monomers.
    monomer_amns = set()
    monomer_mapping = {}
    for sub, identity in subs:
        for atom in sub.GetAtoms():
            amn = atom.GetAtomMapNum()
            if amn > 0:
                monomer_amns.add(amn)
                monomer_mapping[amn] = identity

    # Filter reaction tree for nodes that contain all core type AMNs.
    mols = []
    for _, props in reaction_tree.items():
        mol = Chem.MolFromSmiles(props["smiles"])
        amns = set()
        for atom in mol.GetAtoms():
            amn = atom.GetAtomMapNum()
            if amn > 0: 
                amns.add(amn)
        if monomer_amns.issubset(amns):
            mols.append(mol)

    # Get num of cycles for every filtered out molecule.
    mols_with_num_cycles = []
    for mol in mols:
        ssr = Chem.GetSymmSSSR(mol)
        num_cycles = len(ssr)
        mols_with_num_cycles.append((mol, num_cycles))

    # Get all molecules with the minimum number of cycles.
    min_num_cycles = min([num_cycles for _, num_cycles in mols_with_num_cycles])
    mols = [mol for mol, num_cycles in mols_with_num_cycles if num_cycles == min_num_cycles]

    # If there are no molecules left, return an empty list.
    found_seqs = []

    # Loop over left over molecules and create a monomer graph for each.
    for mol in mols:
        # Create a monomer graph.
        temp_monomer_graph = nx.Graph()
        for atom in mol.GetAtoms():
            amn = atom.GetAtomMapNum()
            if amn > 0 and amn in monomer_amns:
                temp_monomer_graph.add_node(amn)
        for bond in mol.GetBonds():
            begin_amn = bond.GetBeginAtom().GetAtomMapNum()
            end_amn = bond.GetEndAtom().GetAtomMapNum()
            if (
                begin_amn > 0 
                and end_amn > 0
                and begin_amn in monomer_amns
                and end_amn in monomer_amns
            ):
                temp_monomer_graph.add_edge(begin_amn, end_amn)

        # Check if temp_monomer_graph is connected. If not, skip it.
        # This implies that a part of the linear molecule was not identified.
        if not nx.is_connected(temp_monomer_graph):
            logger.debug(f"Found linearized monomer graph that is not connected.")
            continue
        
        contracted_nodes = {}
        for sub, _ in subs:
            sub_amns = []
            for atom in sub.GetAtoms():
                amn = atom.GetAtomMapNum()
                if amn > 0:
                    sub_amns.append(amn)
            
            contracted_nodes[sub_amns[0]] = sub_amns
            # temp_monomer_graph.add_node(amns[0])
            for sub_amn in sub_amns[1:]:
                temp_monomer_graph = nx.contracted_nodes(
                    temp_monomer_graph, 
                    sub_amns[0], 
                    sub_amn, 
                    self_loops=False
                )
        
        # Check if the monomer graph is linear. If not, skip it.
        # The monomer graph is linear if it has one less edge than nodes.
        if temp_monomer_graph.number_of_edges() != temp_monomer_graph.number_of_nodes() - 1:
            logger.debug(f"Found linearized monomer graph that is not linear.")
            continue

        # Check which side is the start and which side is the end.
        # First get both ends, these are the nodes that only have one edge.
        ends = [node for node, degree in temp_monomer_graph.degree() if degree == 1]
        
        # As sanity check, there should be exactly two ends.
        if len(ends) != 2:
            logger.debug(f"Found linearized monomer graph with more than two ends.")
            continue

        # Get the AMNS of both ends.
        amns_a = contracted_nodes[ends[0]]
        amns_b = contracted_nodes[ends[1]]

        # COOH highlights the end of the molecule. Find the end that has a COOH.
        # First find all COOH groups in the linearized molecule.
        coohs = []
        pattern = Chem.MolFromSmarts("*-[C](-[OH])=[O]")
        for match_atom_inds in mol.GetSubstructMatches(pattern):
            cooh = []
            for atom_ind in match_atom_inds:
                atom = mol.GetAtomWithIdx(atom_ind)
                amn = atom.GetAtomMapNum()
                if amn > 0:
                    cooh.append(amn)
            coohs.append(cooh)

        # Check if a and b have a COOH group.
        a_has_cooh = any([set(cooh).issubset(amns_a) for cooh in coohs])
        b_has_cooh = any([set(cooh).issubset(amns_b) for cooh in coohs])

        # If both have a COOH group, skip this molecule.
        if a_has_cooh and b_has_cooh:
            logger.debug(f"Found linearized monomer graph with two ends with COOH.")
            continue

        # If none have a COOH group, skip this molecule.
        if not a_has_cooh and not b_has_cooh:
            logger.debug(f"Found linearized monomer graph with two ends without COOH.")
            continue

        # If a has a COOH group, a is the end.
        if a_has_cooh:
            start_monomer = ends[1]
            end_monomer = ends[0]
        else:
            start_monomer = ends[0]
            end_monomer = ends[1]

        # Get the path between start and end in temp_monomer_graph.
        path = nx.shortest_path(temp_monomer_graph, start_monomer, end_monomer)

        # Get identities of the monomers in the path.
        identities = [monomer_mapping[amn] for amn in path]

        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(f"Found linearized molecule: {Chem.MolToSmiles(mol)}")
            logger.debug(f"Found monomer sequence for molecule: {identities}")

        found_seqs.append(dict(
            smiles=Chem.MolToSmiles(mol),
            motif_code=identities
        ))
    
    return found_seqs
