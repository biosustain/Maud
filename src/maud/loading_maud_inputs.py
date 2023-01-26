"""Provides function load_maud_input."""

import os

import toml

from maud.data_model.experiment import parse_experiment
from maud.data_model.maud_config import MaudConfig
from maud.data_model.maud_init import InitInput
from maud.data_model.maud_input import MaudInput
from maud.data_model.prior_input import PriorInput
from maud.parsing_kinetic_models import parse_kinetic_model


def load_maud_input(data_path: str) -> MaudInput:
    """
    Load an MaudInput object from a data path.

    :param filepath: path to directory containing input toml file
    :param mode: determines which experiments will be included,
    defined in the `experimental_setup` as "sample" and/or "predict".

    """
    # get config
    config_path = os.path.join(data_path, "config.toml")
    config = MaudConfig(**toml.load(config_path))
    # loading
    kinetic_model_path = os.path.join(data_path, config.kinetic_model_file)
    experiments_path = os.path.join(data_path, config.experiments_file)
    raw_kinetic_model = toml.load(kinetic_model_path)
    prior_input_path = os.path.join(data_path, config.priors_file)
    if config.user_inits_file is not None:
        init_input_path = os.path.join(data_path, config.user_inits_file)
        init_input = InitInput(**toml.load(init_input_path))
    else:
        init_input = InitInput()
    # parsing
    kinetic_model = parse_kinetic_model(raw_kinetic_model)
    experiments = [
        parse_experiment(e) for e in toml.load(experiments_path)["experiment"]
    ]
    prior_input = PriorInput(**toml.load(prior_input_path))
    return MaudInput(
        config=config,
        kinetic_model=kinetic_model,
        prior_input=prior_input,
        experiments=experiments,
        init_input=init_input,
    )
