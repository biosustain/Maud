"""Provides dataclass MaudInput."""

from dataclasses import field
from typing import Dict, List, Optional

import pandas as pd
from pydantic.dataclasses import dataclass

from maud.data_model.kinetic_model import KineticModel
from maud.data_model.maud_config import MaudConfig
from maud.data_model.maud_init import InitDict
from maud.data_model.measurement_set import MeasurementSet
from maud.data_model.prior_set import PriorSet, UserPriorInput
from maud.data_model.stan_variable_set import StanVariableSet
from maud.data_model.stan_input import StanInput
from maud.getting_inits import get_inits
from maud.getting_priors import get_prior_set
from maud.getting_stan_variables import get_stan_variable_set
from maud.getting_stan_inputs import get_stan_input


class MIConfig:
    arbitrary_types_allowed = True


@dataclass(config=MIConfig)
class MaudInput:
    """Everything that is needed to run Maud.

    :param kinetic_model: a KineticModel object
    :param priors: a dictionary mapping prior types to lists of Prior objects
    :param stan_coords: a StanCoordSet object
    :param measurement_set: a list of Measurement objects
    :param inits: a dictionary of initial parameter values
    """

    config: MaudConfig
    kinetic_model: KineticModel
    user_priors: UserPriorInput
    measurements: MeasurementSet
    user_inits: Optional[pd.DataFrame]
    stan_variable_set: StanVariableSet = field(init=False)
    priors: PriorSet = field(init=False)
    inits: InitDict = field(init=False)
    stan_input: StanInput = field(init=False)

    def __post_init__(self):
        self.stan_variable_set = get_stan_variable_set(
            self.kinetic_model, self.measurements
        )
        self.priors = get_prior_set(self.user_priors, self.stan_variable_set)
        self.inits = get_inits(self.priors, self.user_inits)
        self.stan_input = get_stan_input(
            self.measurements, self.priors, self.kinetic_model, self.config
        )
