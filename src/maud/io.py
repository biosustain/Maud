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

"""Functions for loading MaudInput objects.

(and at some point in the future also saving MaudOutput objects)

"""

import os

import pandas as pd
import toml

from maud.data_model.maud_input import MaudInput
from maud.data_model.prior_set import UserPriorInput
from maud.parsing import parse_config, parse_kinetic_model, parse_measurement_set
from maud.utils import check_is_df


def load_maud_input(data_path: str) -> MaudInput:
    """
    Load an MaudInput object from a data path.

    :param filepath: path to directory containing input toml file
    :param mode: determines which experiments will be included,
    defined in the `biological_config` as "sample" and/or "predict".

    """
    raw_path = os.path.join(data_path, "config.toml")
    raw = toml.load(raw_path)
    assert isinstance(raw, dict), f"Failed to load input from {raw_path}"
    config = parse_config(raw)
    kinetic_model_path = os.path.join(data_path, config.kinetic_model_file)
    biological_config_path = os.path.join(
        data_path, config.biological_config_file
    )
    raw_kinetic_model = dict(toml.load(kinetic_model_path))
    raw_biological_config = dict(toml.load(biological_config_path))
    measurements_path = os.path.join(data_path, config.measurements_file)
    priors_path = os.path.join(data_path, config.priors_file)
    kinetic_model = parse_kinetic_model(raw_kinetic_model)
    raw_measurements = check_is_df(pd.read_csv(measurements_path))
    main_prior_table = check_is_df(pd.read_csv(priors_path))
    measurement_set = parse_measurement_set(
        raw_measurements, raw_biological_config
    )
    if config.multivariate_dgf_priors:
        dgf_loc = check_is_df(
            pd.read_csv(
                os.path.join(data_path, str(config.dgf_mean_file)),
                index_col="metabolite",
                squeeze=False,
            )
        )["prior_mean_dgf"]
        dgf_cov = check_is_df(
            pd.read_csv(
                os.path.join(data_path, str(config.dgf_covariance_file)),
                index_col="metabolite",
            )
        )
    else:
        dgf_loc = None
        dgf_cov = None
    user_priors = UserPriorInput(main_prior_table, dgf_loc, dgf_cov)
    if config.user_inits_file is not None:
        user_inits_path = os.path.join(data_path, config.user_inits_file)
        user_inits = check_is_df(pd.read_csv(user_inits_path))
    else:
        user_inits = None
    return MaudInput(
        config=config,
        kinetic_model=kinetic_model,
        user_priors=user_priors,
        measurements=measurement_set,
        user_inits=user_inits,
    )
