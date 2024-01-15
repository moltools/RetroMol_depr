"""
Retrieval module.
"""
import typing as ty 

import numpy as np
from sklearn.neighbors import KDTree

def get_nearest_neighbors(fingerprints: np.ndarray, query: np.ndarray, k: int) -> ty.Tuple[np.ndarray, np.ndarray]: 
    """
    Use KD tree to get nearest neighbors.

    :param np.ndarray fingerprints: Fingerprints.
    :param np.ndarray query: Query.
    :param int k: Number of nearest neighbors.
    :return: Nearest neighbors with distances.
    :rtype: ty.Tuple[np.ndarray, np.ndarray]
    """
    tree = KDTree(fingerprints, leaf_size=2)
    dist, ind = tree.query(query, k=k)
    return dist[0], ind[0]