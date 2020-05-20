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

from maud.data_model import KineticModel


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


def codify(lx: Iterable[str]) -> Dict[str, int]:
    """Turn an iterable of strings into a dictionary mapping them to integer indexes."""
    return dict(zip(lx, range(1, len(lx) + 1)))


def get_metabolite_codes(kinetic_model: KineticModel) -> Dict[str, int]:
    """Get a dictionary mapping metabolite ids to integer indexes.

    :param kinetic_model: A KineticModel object
    """

    return codify(kinetic_model.metabolites.keys())


def get_enzyme_codes(kinetic_model: KineticModel) -> Dict[str, int]:
    """Get a dictionary mapping enzyme ids to integer indexes.

    :param kinetic_model: A KineticModel object
    """

    enzyme_ids = []
    for _, rxn in kinetic_model.reactions.items():
        for enz_id, _ in rxn.enzymes.items():
            enzyme_ids.append(enz_id)
    return codify(enzyme_ids)


def get_kinetic_parameter_codes(kinetic_model: KineticModel) -> Dict[str, int]:
    """Get a dictionary mapping kinetic_parameter ids to integer indexes.

    :param kinetic_model: A KineticModel object
    """

    parameter_ids = []
    for _, rxn in kinetic_model.reactions.items():
        for enz_id, enz in rxn.enzymes.items():
            for par_id, par in enz.parameters.items():
                if not par.is_thermodynamic:
                    parameter_ids.append(enz_id + "_" + par_id)
    return codify(parameter_ids)


def get_thermo_codes(kinetic_model: KineticModel) -> Dict[str, int]:
    """Get a dictionary mapping enzyme ids to integer indexes.

    :param kinetic_model: A KineticModel object
    """
    parameter_ids = []
    for _, rxn in kinetic_model.reactions.items():
        for enz_id, enz in rxn.enzymes.items():
            for par_id, par in enz.parameters.items():
                if par.is_thermodynamic:
                    parameter_ids.append(enz_id + "_" + par_id)
    return codify(parameter_ids)
