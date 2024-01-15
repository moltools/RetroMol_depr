"""
Module for biosynthetic fingerprinting.
"""
import math 
import typing as ty
from enum import Enum

import numpy as np

class CompoundClassMapping(Enum):
    """
    Enum class that maps motif class to coordinates on unit circle.

    For more information on the compound classes, see:
    https://srobinson.shinyapps.io/AdenylPred/#section-substrate-groups 
    """
    A =                     ( 1.000,  0.000)
    B =                     ( 0.841,  0.541)
    C =                     ( 0.415,  0.910)
    D =                     (-0.142,  0.990)
    PolarAndCharged =       (-0.655,  0.756)
    SmallHydrophobic =      (-0.959,  0.282)
    SmallNonHydrophobic =   (-0.655, -0.756)
    Tiny =                  ( 0.841, -0.541)
    Bulky =                 (-0.959, -0.282)
    CyclicAliphatic =       (-0.142, -0.990)
    SulfurContaining =      ( 0.415, -0.910)
    Undefined =              ( 0.000,  0.000)

    def get_mapping(item: str) -> ty.Tuple[float, float]:
        """
        Get mapping on unit circle for motif class.
        
        :param str item: Item.
        :return: Mapping.
        :rtype: ty.Tuple[float, float]
        :raises ValueError: If item is not a valid motif class.
        """
        return CompoundClassMapping[item].value
    
def get_biosynthetic_fingerprint(seq: ty.List[str]) -> np.ndarray:
    """
    Create fingerprint from sequence, based on Chaos Game Representation.

    :param ty.List[str] seq: List of motif classes.
    :return: Biosynthetic fingerprint (N=3600).
    :rtype: np.ndarray
    """
    # Make empty 3600 count array.
    fingerprint = np.zeros(3600, dtype=np.int8)

    def find_bin(x: float, y: float) -> int:
        """
        Find bin index for given coordinates.

        :param float x: x-coordinate.
        :param float y: y-coordinate.
        :return: Bin index.
        :rtype: int
        """
        # Get clockwise angle from x-axis.
        angle = math.degrees(math.atan2(y, x))
        if angle < 0: angle += 360
        angle = math.ceil(angle)

        # Get distance to origin and multiply by 10 do get distance to 
        # origin as integer between 0 and 10.
        dist = math.ceil(math.sqrt(x**2 + y**2) * 10)

        # Get bin index.
        bin_index = angle * dist

        return bin_index
    
    prev_loc = (0, 0)
    for motif in seq:
        next_loc = CompoundClassMapping.get_mapping(motif)        
        middle_loc = ((prev_loc[0] + next_loc[0]) / 2, (prev_loc[1] + next_loc[1]) / 2)
        bin_index = find_bin(*middle_loc)
        fingerprint[bin_index] += 1
        prev_loc = middle_loc

    # Normalize fingerprint.
    fingerprint = fingerprint / np.linalg.norm(fingerprint)

    return fingerprint