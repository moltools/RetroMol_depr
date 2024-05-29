# -*- coding: utf-8 -*-

"""This module contains functions for sequence alignment."""

import re
import typing as ty

from versalign.motif import Motif
from versalign.sequence import Sequence


class PolyketideMotif(Motif):
    """A polyketide motif."""

    def __init__(self, type: ty.Optional[str] = None, decoration: ty.Optional[int] = None) -> None:
        """Initialize a polyketide motif.

        :param type: Type of polyketide.
        :type type: ty.Optional[str]
        :param decoration: Decoration of polyketide.
        :type decoration: ty.Optional[int]
        """
        super().__init__()

        self.type = type
        self.decoration = decoration

    def __eq__(self, other: ty.Any) -> bool:
        """Check if two polyketide motifs are equal.

        :param other: Other polyketide motif.
        :type other: ty.Any
        :return: True if equal, False otherwise.
        :rtype: bool
        """
        if isinstance(other, PolyketideMotif):
            return self.type == other.type and self.decoration == other.decoration

        return False

    def __str__(self) -> str:
        """Return string representation of polyketide motif.

        :return: String representation of polyketide motif.
        :rtype: str
        """
        if self.type is not None and self.decoration is not None:
            return f"{self.type}{self.decoration}"

        elif self.type is not None:
            return f"{self.type}"

        return "?"


class PeptideMotif(Motif):
    """A peptide motif."""

    def __init__(self, source: ty.Optional[str] = None, cid: ty.Optional[str] = None) -> None:
        """Initialize a peptide motif.

        :param type: Type of peptide.
        :type type: ty.Optional[str]
        :param sequence: Sequence of peptide.
        :type sequence: ty.Optional[str]
        """
        super().__init__()

        self.source = source
        self.cid = cid

    def __eq__(self, other: ty.Any) -> bool:
        """Check if two peptide motifs are equal.

        :param other: Other peptide motif.
        :type other: ty.Any
        :return: True if equal, False otherwise.
        :rtype: bool
        """
        if isinstance(other, PeptideMotif):
            return self.source == other.source and self.cid == other.cid

        return False

    def __str__(self) -> str:
        """Return string representation of peptide motif.

        :return: String representation of peptide motif.
        :rtype: str
        """
        if self.source is not None and self.cid is not None:
            return f"{self.source}:{self.cid}"

        return "?"


def sequence_from_motif_string_list(name: str, motif_string_list: ty.List[str]) -> Sequence:
    """Create a sequence from a list of motif strings.

    :param name: Name of sequence.
    :type name: str
    :param string_list: List of motif strings.
    :type string_list: ty.List[str]
    :return: Sequence.
    :rtype: Sequence
    """
    motifs = []

    for motif_string in motif_string_list:
        if match := re.match(r"polyketide\|([A-D])(\d{1,2})", motif_string):
            motif_type = "polyketide"
            calculated = True
            polyketide_type = match.group(1)
            polyketide_decoration_type = int(match.group(2))

            motif = PolyketideMotif(type=polyketide_type, decoration=polyketide_decoration_type)
            motifs.append(motif)

        elif match := re.match(r"peptide\|(\w+)\|(.+)", motif_string):
            motif_type = "peptide"
            calculated = True
            peptide_source = match.group(1)
            peptide_cid = match.group(2)

            motif = PeptideMotif(source=peptide_source, cid=peptide_cid)
            motifs.append(motif)

        else:
            raise ValueError(f"Unknown motif: {motif_string}")

    return Sequence(name, motifs)
