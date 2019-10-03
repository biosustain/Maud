# Copyright (C) 2019 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
