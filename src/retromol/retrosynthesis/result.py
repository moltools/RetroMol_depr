# -*- coding: utf-8 -*-

"""This module contains the dataclass for storing the results of a RetroMol run."""

import json
import typing as ty
from dataclasses import dataclass

from rdkit import Chem


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
