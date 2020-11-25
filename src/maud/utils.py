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

from typing import Dict, Iterable, List

import numpy as np
from scipy.stats import norm

from maud.data_model import KineticModel


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


def get_lognormal_parameters_from_quantiles(x1, p1, x2, p2):
    """Find parameters for a lognormal distribution from two quantiles.

    i.e. get mu and sigma such that if X ~ lognormal(mu, sigma), then pr(X <
    x1) = p1 and pr(X < x2) = p2.

    """
    logx1 = np.log(x1)
    logx2 = np.log(x2)
    denom = norm.ppf(p2) - norm.ppf(p1)
    sigma = (logx2 - logx1) / denom
    mu = (logx1 * norm.ppf(p2) - logx2 * norm.ppf(p1)) / denom
    return mu, sigma


def get_normal_parameters_from_quantiles(x1, p1, x2, p2):
    """Find parameters for a normal distribution from two quantiles.

    i.e. get mu and sigma such that if X ~ normal(mu, sigma), then pr(X <
    x1) = p1 and pr(X < x2) = p2.

    """
    denom = norm.ppf(p2) - norm.ppf(p1)
    sigma = (x2 - x1) / denom
    mu = (x1 * norm.ppf(p2) - x2 * norm.ppf(p1)) / denom
    return mu, sigma


def standardise_list(l: List[float], mean: float, scale: float) -> List[float]:
    """Subtract the mean and divide by the sd."""
    return ((np.array(l) - np.array(mean)) / np.array(scale)).tolist()
