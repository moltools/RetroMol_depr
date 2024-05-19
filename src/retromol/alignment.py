# -*- coding: utf-8 -*-

"""This module contains functions for sequence alignment."""

import re
import typing as ty

from versalign.motif import Motif 
from versalign.sequence import Sequence
from versalign.pairwise import PairwiseAlignment, align_pairwise


class PolyketideMotif(Motif):
    """A polyketide motif."""

    def __init__(
        self, 
        type: ty.Optional[str] = None,
        decoration: ty.Optional[int] = None
    ) -> None:
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

    def __init__(
        self, 
        source: ty.Optional[str] = None, 
        cid: ty.Optional[str] = None
    ) -> None:
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
