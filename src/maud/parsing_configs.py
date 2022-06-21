"""Functions for turning raw data into Maud objects."""

from typing import Dict

from maud.data_model.maud_config import MaudConfig, ODEConfig


def parse_config(raw: Dict) -> MaudConfig:
    """Get a MaudConfig object from the result of toml.load.

    :param raw: result of running toml.load on a suitable file
    """
    init_kwargs = raw.copy()
    if "ode_config" in init_kwargs.keys():
        init_kwargs["ode_config"] = ODEConfig(**init_kwargs["ode_config"])
    return MaudConfig(**init_kwargs)
