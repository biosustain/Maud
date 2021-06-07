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
import pandas as pd
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




def get_input_template(km, raw_measurements):
    """Extact parameters from an infd object."""

    from maud.io import get_stan_coords

    class Input_Coords:
        """Defines parameters with associated coordinate sets.

        :param id: id of the parameter.
        :param coords: dictionary.
        """

        def __init__(self, id: str, coords: Dict[str, List[str]]):
            self.id = id
            self.coords = coords

    def get_2d_coords(coords_1, coords_2):
        """Return unpacked coordinates for 2-D indexing."""
        set_of_coords = []
        for c1 in coords_1:
            for c2 in coords_2:
                set_of_coords.append((c1, c2))
        return list(zip(*set_of_coords)) if len(set_of_coords) > 0 else ([], [])

    scs = get_stan_coords(km, raw_measurements)

    list_of_input_inits = [
        Input_Coords(
            id="km",
            coords={"enzyme_id": scs.km_enzs, "mic_id": scs.km_mics},
        ),
        Input_Coords(
            id="drain",
            coords={
                "drain_id": get_2d_coords(scs.drains, scs.experiments)[0],
                "experiment_id": get_2d_coords(scs.drains, scs.experiments)[1],
            },
        ),
        Input_Coords(id="ki", coords={"enzyme_id": scs.ci_enzs, "mic_id": scs.ci_mics}),
        Input_Coords(
            id="diss_t", coords={"enzyme_id": scs.ai_enzs, "mic_id": scs.ai_mics}
        ),
        Input_Coords(
            id="diss_r", coords={"enzyme_id": scs.aa_enzs, "mic_id": scs.aa_mics}
        ),
        Input_Coords(
            id="transfer_constant", coords={"enzyme_id": scs.allosteric_enzymes}
        ),
        Input_Coords(id="kcat", coords={"enzyme_id": scs.enzymes}),
        Input_Coords(id="kcat_phos", coords={"phos_enz_id": scs.phos_enzs}),
        Input_Coords(
            id="conc_unbalanced",
            coords={
                "mic_id": get_2d_coords(scs.unbalanced_mics, scs.experiments)[0],
                "experiment_id": get_2d_coords(scs.unbalanced_mics, scs.experiments)[1],
            },
        ),
        Input_Coords(
            id="conc_enzyme",
            coords={
                "enzyme_id": get_2d_coords(scs.enzymes, scs.experiments)[0],
                "experiment_id": get_2d_coords(scs.enzymes, scs.experiments)[1],
            },
        ),
        Input_Coords(
            id="conc_phos",
            coords={
                "phos_enz_id": get_2d_coords(scs.phos_enzs, scs.experiments)[0],
                "experiment_id": get_2d_coords(scs.phos_enzs, scs.experiments)[1],
            },
        ),
        Input_Coords(id="dgf", coords={"metabolite_id": scs.metabolites}),
    ]

    init_dataframe = pd.DataFrame(
        columns=[
            "parameter_name",
            "experiment_id",
            "metabolite_id",
            "mic_id",
            "enzyme_id",
            "phos_enz_id",
            "drain_id",
            "location",
            "scale",
            "pct1",
            "pct99",
        ]
    )

    for par in list_of_input_inits:
        par_dataframe = pd.DataFrame.from_dict(par.coords)
        par_dataframe["parameter_name"] = par.id
        init_dataframe = init_dataframe.append(par_dataframe, ignore_index=True)

    return init_dataframe
