"""General purpose utility functions."""

import os
from typing import Dict, Iterable

import numpy as np


def match_string_to_file(s: str, path: str) -> bool:
    """Check if a string is the same as the contents of a file."""
    if os.path.exists(path):
        return open(path, "r").read() == s
    else:
        return False


def sem_pct_to_lognormal_sigma(sem_pct, mean, n=3):
    """Get the lognormal sigma parameter for an sem percentage."""
    sem = sem_pct / 100 * mean
    s = sem * np.sqrt(n)
    return np.sqrt(np.log(1 + (s ** 2) / (mean ** 2)))


def codify(l: Iterable[str]) -> Dict[str, int]:
    """Turn an iterable of strings into a dictionary mapping them to integer indexes."""
    return dict(zip(l, range(1, len(l) + 1)))
