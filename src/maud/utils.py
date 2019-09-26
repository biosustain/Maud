# Copyright (c) 2019, Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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
