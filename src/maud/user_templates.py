# Copyright (C) 2021 Novo Nordisk Foundation Center for Biosustainability,
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

"""Functions and classes to generate prior template stored here."""

from typing import Dict, List, Optional

import pandas as pd

from maud.analysis import join_list_of_strings
from maud.io import get_stan_coords


PRIOR_FILE_COLUMNS = [
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

INIT_FILE_COLUMNS = [
    "parameter_name",
    "experiment_id",
    "metabolite_id",
    "mic_id",
    "enzyme_id",
    "phos_enz_id",
    "drain_id",
    "value",
]


class Input_Coords:
    """Defines parameters with associated coordinate sets.

    :param id: id of the parameter.
    :param coords: dictionary of user input keys with corresponding list of ordered ids.
    :param infd_coord_list: ordered list that matches the coords to the infd object.
    :param linking_list: matches composite strings to the infd coords.
    """

    def __init__(
        self,
        id: str,
        coords: Dict[str, List[str]],
        infd_coord_list: List[str],
        linking_list: Optional[Dict[str, List[str]]] = None,
    ):
        self.id = id
        self.coords = coords
        self.infd_coord_list = infd_coord_list
        self.linking_list = linking_list


def get_2d_coords(coords_1, coords_2):
    """Return unpacked coordinates for 2-D indexing."""
    set_of_coords = []
    for c1 in coords_1:
        for c2 in coords_2:
            set_of_coords.append((c1, c2))
    return list(zip(*set_of_coords)) if len(set_of_coords) > 0 else ([], [])


def get_parameter_coords(scs):
    """Define parameter coordinates for stan and infd objects."""
    return [
        Input_Coords(
            id="km",
            coords={"enzyme_id": scs.km_enzs, "mic_id": scs.km_mics},
            infd_coord_list=["kms"],
            linking_list={"kms": join_list_of_strings(scs.km_enzs, scs.km_mics)},
        ),
        Input_Coords(
            id="drain",
            coords={
                "drain_id": get_2d_coords(scs.drains, scs.experiments)[0],
                "experiment_id": get_2d_coords(scs.drains, scs.experiments)[1],
            },
            infd_coord_list=["drains", "experiments"],
        ),
        Input_Coords(
            id="ki",
            coords={"enzyme_id": scs.ci_enzs, "mic_id": scs.ci_mics},
            infd_coord_list=["kis"],
            linking_list={"kis": join_list_of_strings(scs.ci_enzs, scs.ci_mics)},
        ),
        Input_Coords(
            id="diss_t",
            coords={"enzyme_id": scs.ai_enzs, "mic_id": scs.ai_mics},
            infd_coord_list=["diss_ts"],
            linking_list={"diss_ts": join_list_of_strings(scs.ai_enzs, scs.ai_mics)},
        ),
        Input_Coords(
            id="diss_r",
            coords={"enzyme_id": scs.aa_enzs, "mic_id": scs.aa_mics},
            infd_coord_list=["diss_rs"],
            linking_list={"diss_rs": join_list_of_strings(scs.aa_enzs, scs.aa_mics)},
        ),
        Input_Coords(
            id="transfer_constant",
            coords={"enzyme_id": scs.allosteric_enzymes},
            infd_coord_list=["allosteric_enzymes"],
        ),
        Input_Coords(
            id="kcat",
            coords={"enzyme_id": scs.enzymes},
            infd_coord_list=["enzymes"],
        ),
        Input_Coords(
            id="kcat_phos",
            coords={"phos_enz_id": scs.phos_enzs},
            infd_coord_list=["phos_enzs"],
        ),
        Input_Coords(
            id="conc_unbalanced",
            coords={
                "mic_id": get_2d_coords(scs.unbalanced_mics, scs.experiments)[0],
                "experiment_id": get_2d_coords(scs.unbalanced_mics, scs.experiments)[1],
            },
            infd_coord_list=["unbalanced_mics", "experiments"],
        ),
        Input_Coords(
            id="conc_enzyme",
            coords={
                "enzyme_id": get_2d_coords(scs.enzymes, scs.experiments)[0],
                "experiment_id": get_2d_coords(scs.enzymes, scs.experiments)[1],
            },
            infd_coord_list=["enzymes", "experiments"],
        ),
        Input_Coords(
            id="conc_phos",
            coords={
                "phos_enz_id": get_2d_coords(scs.phos_enzs, scs.experiments)[0],
                "experiment_id": get_2d_coords(scs.phos_enzs, scs.experiments)[1],
            },
            infd_coord_list=["phos_enzs", "experiments"],
        ),
        Input_Coords(
            id="dgf",
            coords={"metabolite_id": scs.metabolites},
            infd_coord_list=["metabolites"],
        ),
        Input_Coords(id="keq", coords={"edges": scs.edges}, infd_coord_list=["edges"]),
    ]


def get_prior_template(km, raw_measurements, experiments, mode):
    """Get prior dataframe from KineticModel and Measurements."""

    scs = get_stan_coords(km, raw_measurements, experiments, mode)
    list_of_input_inits = get_parameter_coords(scs)
    prior_dataframe = pd.DataFrame(columns=PRIOR_FILE_COLUMNS)
    for par in list_of_input_inits:
        par_dataframe = pd.DataFrame.from_dict(par.coords)
        par_dataframe["parameter_name"] = par.id
        prior_dataframe = prior_dataframe.append(par_dataframe, ignore_index=True)

    return prior_dataframe


def get_inits_from_draw(infd, mi, chain, draw, warmup):
    """Extact parameters from an infd object."""

    scs = mi.stan_coords
    list_of_input_inits = get_parameter_coords(scs)
    init_dataframe = pd.DataFrame(columns=INIT_FILE_COLUMNS)
    if warmup == 1:
        infd_parameters = list(infd.warmup_posterior.variables.keys())
    else:
        infd_parameters = list(infd.posterior.variables.keys())
    for par in list_of_input_inits:
        if par.id in infd_parameters:
            if warmup == 1:
                value_dataframe = (
                    infd.warmup_posterior[par.id][chain][draw]
                    .to_dataframe()
                    .reset_index()
                )

            else:
                value_dataframe = (
                    infd.posterior[par.id][chain][draw].to_dataframe().reset_index()
                )
            par_dataframe = pd.DataFrame.from_dict(par.coords)
            if par.linking_list:
                par_dataframe["linking_list"] = list(par.linking_list.values())[0]
                par_dataframe = par_dataframe.merge(
                    value_dataframe,
                    left_on="linking_list",
                    right_on=par.infd_coord_list,
                )
            else:
                par_dataframe = par_dataframe.merge(
                    value_dataframe,
                    left_on=list(par.coords.keys()),
                    right_on=par.infd_coord_list,
                )
            par_dataframe["parameter_name"] = par.id
            par_dataframe["value"] = par_dataframe[par.id]
            init_column_list = ["parameter_name"] + list(par.coords.keys()) + ["value"]
            init_dataframe = init_dataframe.append(par_dataframe[init_column_list])
    return init_dataframe
