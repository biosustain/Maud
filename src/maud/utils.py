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

from typing import Dict, Iterable

import numpy as np
import sympy as sp
from scipy.stats import norm


def codify(lx: Iterable[str]) -> Dict[str, int]:
    """Turn an iterable of strings into a dictionary mapping them to integer indexes."""
    return dict(zip(lx, range(1, len(lx) + 1)))


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


def get_null_space(a, rtol=1e-5):
    """Calulate the null space of a matrix."""
    u, s, v = np.linalg.svd(a)
    rank = (s > rtol * s[0]).sum()
    return v[rank:].T.copy()


def get_rref(mat):
    """Return reduced row echelon form of a matrix."""
    return sp.Matrix(mat).rref(iszerofunc=lambda x: abs(x) < 1e-10)[0]


def export_scaled_params_from_draws(infd, chain, draw):
    """Extact log centered parameters from an infd object."""
    list_of_input_inits = [
        "log_km_z",
        "drain_z",
        "log_ki_z",
        "log_dissociation_constant_t_z",
        "log_dissociation_constant_r_z",
        "log_transfer_constant_z",
        "log_kcat_z",
        "log_phos_kcat_z",
        "log_conc_unbalanced_z",
        "log_enzyme_z",
        "log_phos_conc_z",
        "formation_energy_z",
    ]

    input_dict = {
        par_name: infd.posterior[par_name][chain, draw].values.tolist()
        for par_name in list_of_input_inits
        if par_name in infd.posterior.variables.keys()
    }
    return input_dict


def export_params_from_draw(infd, chain, draw):
    """Extact parameters from an infd object."""
    list_of_input_inits = [
        "km",
        "drain",
        "ki",
        "dissociation_constant_t",
        "dissociation_constant_r",
        "transfer_constant",
        "kcat",
        "phos_enzyme_kcat",
        "conc_unbalanced",
        "enzyme",
        "phos_enzyme_conc",
        "formation_energy",
    ]

    input_dict = {
        par_name: infd.posterior[par_name][chain, draw].values.tolist()
        for par_name in list_of_input_inits
        if par_name in infd.posterior.variables.keys()
    }
    return input_dict
