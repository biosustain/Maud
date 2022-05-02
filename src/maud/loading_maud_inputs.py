"""Provides function load_maud_input."""

import os

import toml

from maud.data_model.maud_input import MaudInput
from maud.data_model.prior_set import UserPriorInput
from maud.parsing_configs import parse_config
from maud.parsing_kinetic_models import parse_kinetic_model
from maud.parsing_measurements import parse_measurement_set
from maud.utils import load_df


def load_maud_input(data_path: str) -> MaudInput:
    """
    Load an MaudInput object from a data path.

    :param filepath: path to directory containing input toml file
    :param mode: determines which experiments will be included,
    defined in the `experimental_setup` as "sample" and/or "predict".

    """
    # get config
    config_path = os.path.join(data_path, "config.toml")
    config = parse_config(toml.load(config_path))
    # loading
    kinetic_model_path = os.path.join(data_path, config.kinetic_model_file)
    bio_config_path = os.path.join(data_path, config.experimental_setup_file)
    raw_kinetic_model = toml.load(kinetic_model_path)
    raw_bio_config = toml.load(bio_config_path)
    measurements_path = os.path.join(data_path, config.measurements_file)
    priors_path = os.path.join(data_path, config.priors_file)
    raw_measurements = load_df(measurements_path)
    main_prior_table = load_df(priors_path)
    if (
        config.dgf_mean_file is not None
        and config.dgf_covariance_file is not None
    ):
        dgf_mean_path = os.path.join(data_path, config.dgf_mean_file)
        dgf_cov_path = os.path.join(data_path, config.dgf_covariance_file)
        dgf_loc = load_df(
            dgf_mean_path, index_col="metabolite", squeeze=False
        )["prior_mean_dgf"]
        dgf_cov = load_df(dgf_cov_path, index_col="metabolite")
    else:
        dgf_loc = None
        dgf_cov = None
    if config.user_inits_file is not None:
        user_inits_path = os.path.join(data_path, config.user_inits_file)
        user_inits = load_df(user_inits_path)
    else:
        user_inits = None
    # parsing
    kinetic_model = parse_kinetic_model(raw_kinetic_model)
    measurement_set = parse_measurement_set(raw_measurements, raw_bio_config)
    user_priors = UserPriorInput(main_prior_table, dgf_loc, dgf_cov)
    return MaudInput(
        config=config,
        kinetic_model=kinetic_model,
        user_priors=user_priors,
        measurements=measurement_set,
        user_inits=user_inits,
    )
