"""Provides dataclass MaudInput containing Maud needs to run."""

from typing import Optional

import pandas as pd
from pydantic import Field
from pydantic.dataclasses import dataclass

from maud.data_model.kinetic_model import KineticModel
from maud.data_model.maud_config import MaudConfig
from maud.data_model.maud_init import InitDict
from maud.data_model.measurement_set import MeasurementSet
from maud.data_model.prior_set import PriorSet, UserPriorInput
from maud.data_model.stan_input import StanInputTest, StanInputTrain
from maud.data_model.stan_variable_set import StanVariableSet
from maud.getting_inits import get_inits
from maud.getting_priors import get_prior_set
from maud.getting_stan_inputs import get_stan_inputs
from maud.getting_stan_variables import get_stan_variable_set


class MIConfig:
    """Config for MaudInput, allowing it to contain pandas objects."""

    arbitrary_types_allowed = True


@dataclass(config=MIConfig)
class MaudInput:
    """Everything that is needed to run Maud."""

    config: MaudConfig
    kinetic_model: KineticModel
    user_priors: UserPriorInput
    measurements: MeasurementSet
    user_inits: Optional[pd.DataFrame]
    stan_variable_set: StanVariableSet = Field(init=False, exclude=True)
    priors: PriorSet = Field(init=False, exclude=True)
    inits: InitDict = Field(init=False, exclude=True)
    stan_input_train: StanInputTrain = Field(init=False, exclude=True)
    stan_input_test: Optional[StanInputTest] = Field(init=False, exclude=True)

    def __post_init__(self):
        """Add attributes that depend on other ones."""
        self.stan_variable_set = get_stan_variable_set(
            self.kinetic_model, self.measurements
        )
        self.priors = get_prior_set(self.user_priors, self.stan_variable_set)
        self.inits = get_inits(self.priors, self.user_inits, self.measurements)
        self.stan_input_train, self.stan_input_test = get_stan_inputs(
            self.measurements, self.priors, self.kinetic_model, self.config
        )
