# -*- coding: utf-8 -*-

"""This module contains the API for RetroMol."""

import json
import logging
import typing as ty
from dataclasses import dataclass

from rdkit import Chem

from retromol.chem import MolecularPattern, Molecule, ReactionRule
from retromol.graph import reaction_tree_to_digraph, reaction_tree_to_monomer_graph
from retromol.sequencing import parse_modular_natural_product


@dataclass
class Result:
    """Dataclass for storing the results of a RetroMol run.

    :param identifier: The identifier of the molecule.
    :type identifier: str
    :param mol: The molecule.
    :type mol: Chem.Mol
    :param success: Whether the RetroMol run was successful.
    :type success: bool
    :param reaction_tree: The reaction tree.
    :type reaction_tree: ty.Dict[int, ty.List[int]]
    :param applied_reactions: The applied reactions.
    :type applied_reactions: ty.List[str]
    :param monomer_graph: The monomer graph.
    :type monomer_graph: ty.Dict[str, ty.Any]
    :param sequences: The modular natural product sequences.
    :type sequences: ty.List[ty.List[str]]
    """

    identifier: str
    mol: Chem.Mol
    success: bool
    reaction_tree: ty.Dict[int, ty.List[int]] = None
    applied_reactions: ty.List[str] = None
    monomer_graph: ty.Dict[str, ty.Any] = None
    sequences: ty.List[ty.List[str]] = None

    def to_json(self, indent: int = 4) -> str:
        """Convert the result to JSON.

        :param indent: The number of spaces to indent the JSON.
        :type indent: int
        :return: The result as JSON.
        :rtype: str
        """
        return json.dumps(
            {
                "identifier": self.identifier,
                "smiles": Chem.MolToSmiles(self.mol),
                "success": "true" if self.success else "false",
                "reaction_tree": self.reaction_tree,
                "applied_reactions": self.applied_reactions,
                "monomer_graph": self.monomer_graph,
                "sequences": self.sequences,
            },
            indent=indent,
        )

    @classmethod
    def from_json(cls, path: str) -> "Result":
        """Create a Result object from a JSON file.

        :param path: The path to the JSON file.
        :type path: str
        :return: The Result object.
        :rtype: Result
        """
        with open(path, "r", encoding="utf-8") as fo:
            data = json.load(fo)

        reaction_tree = data["reaction_tree"]
        if reaction_tree is not None:
            reaction_tree = {int(k): v for k, v in reaction_tree.items()}

        monomer_graph = data["monomer_graph"]
        if monomer_graph is not None:
            monomer_graph = {int(k): v for k, v in monomer_graph.items()}

        return Result(
            identifier=data["identifier"],
            mol=Chem.MolFromSmiles(data["smiles"]),
            success=True if data["success"] == "true" else False,
            reaction_tree=reaction_tree,
            applied_reactions=data["applied_reactions"],
            monomer_graph=monomer_graph,
            sequences=data["sequences"],
        )

    def has_identified_monomers(self) -> bool:
        """Check if any monomers have been identified.

        :return: Whether any monomers have been identified.
        :rtype: bool
        """
        return any([x["identity"] is not None for _, x in self.monomer_graph.items()])


def parse_reaction_rules(data: ty.Dict[str, str]) -> ty.List[ReactionRule]:
    """Parse reaction rules from JSON.

    :param data: The JSON source.
    :type data: str
    :return: The reaction rules.
    :rtype: ty.List[ReactionRule]
    :raises Exception: If 'name' and/or 'pattern' is not provided for JSON
        item.
    """
    reaction_rules = []
    for item in data:
        name = item.get("name", None)
        pattern = item.get("pattern", None)

        if name is None or pattern is None:
            raise Exception("Reaction item invalid: 'name' and/or 'pattern' is not " "provided.")

        reaction_rule = ReactionRule(name, pattern)
        reaction_rules.append(reaction_rule)

    return reaction_rules


def parse_molecular_patterns(data: str) -> ty.List[MolecularPattern]:
    """Parse molecular patterns from JSON.

    :param data: The JSON source.
    :type data: str
    :return: The molecular patterns.
    :rtype: ty.List[MolecularPattern]
    :raises Exception: If 'name' and/or 'pattern' is not provided for JSON
        item.
    """
    molecular_patterns = []
    for item in data:
        name = item.get("name", None)
        pattern = item.get("pattern", None)

        if name is None or pattern is None:
            raise Exception("Monomer item invalid: 'name' and/or 'pattern' is not " "provided.")

        molecular_pattern = MolecularPattern(name, pattern)
        molecular_patterns.append(molecular_pattern)

    return molecular_patterns


def parse_mol(
    mol: Molecule,
    reactions: ty.List[ReactionRule],
    monomers: ty.List[MolecularPattern],
    neutralize: bool = True,
) -> Result:
    """Parse a molecule with RetroMol.

    :param mol: The molecule.
    :type mol: Molecule
    :param reactions: The reaction rules.
    :type reactions: ty.List[ReactionRule]
    :param monomers: The molecular patterns.
    :type monomers: ty.List[MolecularPattern]
    :param neutralize: Whether to neutralize the molecule.
    :type neutralize: bool
    :return: The result.
    :rtype: Result
    """
    logger = logging.getLogger(__name__)

    name = mol.name

    if neutralize:
        mol.neutralize()

    logger.debug(f"Starting parsing for {name} ...")
    reactant, reaction_tree, reaction_mapping = mol.apply_rules(reactions)
    logger.debug(f"Applied reactions to {name}.")

    applied_reactions = list(
        set([reaction for _, reactions in reaction_tree.items() for reaction in reactions])
    )
    logger.debug(f"Applied reactions: {applied_reactions}")

    logger.debug(f"Starting monomer graph generation for {name} ...")
    reaction_tree = reaction_tree_to_digraph(reaction_tree)
    logger.debug("... Created digraph from reaction tree.")
    monomer_graph, monomer_mapping = reaction_tree_to_monomer_graph(
        mol, reaction_tree, reaction_mapping, monomers
    )
    logger.debug("... Generated monomer graph from digraph.")
    logger.debug(f"Generated monomer graph for {name}.")

    # Reformat reaction_tree and reaction_mapping into one graph.
    logger.debug(f"Reformatting reaction tree for {name} ...")
    new_reaction_tree = {}
    for parent in reaction_tree.nodes:
        children = list(reaction_tree.successors(parent))

        mol = reaction_mapping[parent]
        for atom in mol.GetAtoms():
            amn = atom.GetIsotope()
            atom.SetIsotope(0)
            atom.SetAtomMapNum(amn)
        smiles = Chem.MolToSmiles(mol)

        new_reaction_tree[parent] = {"smiles": smiles, "children": children}

    logger.debug(f"Reformatted reaction tree for {name}.")

    # Reformat monomer_graph and monomer_mapping into one graph.
    logger.debug(f"Reformatting monomer graph for {name} ...")
    new_monomer_graph = {}

    for reaction_tree_node, (monomer_graph_node, identity) in monomer_mapping.items():
        new_monomer_graph[monomer_graph_node] = {
            "reaction_tree_id": reaction_tree_node,
            "identity": identity,
            "neighbors": [],
        }

    for node in monomer_graph:
        neighbors = list(monomer_graph[node])

        if node in new_monomer_graph:
            new_monomer_graph[node]["neighbors"] = neighbors
        else:
            new_monomer_graph[node] = {
                "reaction_tree_id": None,
                "identity": None,
                "neighbors": neighbors,
            }

    logger.debug(f"Reformatted monomer graph for {name}.")

    if logger.isEnabledFor(logging.DEBUG):
        logger.debug(f" Parsing was successful for {name}.")

        for node, items in new_monomer_graph.items():
            if items["identity"] is not None:
                logger.debug(f"Found identity for {node}: {items['identity']}")

    try:
        seqs = parse_modular_natural_product(new_reaction_tree, new_monomer_graph)
    except Exception as e:
        logger.error(f" Error while parsing modular natural product: {e}")
    finally:
        if logger.isEnabledFor(logging.DEBUG):
            if not len(seqs):
                logger.debug("No modular natural product sequence found.")
            for seq in seqs:
                logger.debug(f"Modular natural product sequence: {seq}")

    # Create result object.
    result = Result(
        identifier=name,
        mol=reactant,
        success=True,
        reaction_tree=new_reaction_tree,
        applied_reactions=applied_reactions,
        monomer_graph=new_monomer_graph,
        sequences=seqs,
    )

    return result
