"""Provides dataclass MaudInput containing Maud needs to run."""

from dataclasses import fields
from typing import Dict, List

from pydantic import Field
from pydantic.dataclasses import dataclass

from maud.data_model.experiment import Experiment
from maud.data_model.kinetic_model import KineticModel
from maud.data_model.maud_config import MaudConfig
from maud.data_model.maud_init import InitInput
from maud.data_model.maud_parameter import ParameterSet
from maud.data_model.prior_input import PriorInput
from maud.getting_parameters import get_maud_parameters
from maud.getting_stan_inputs import get_stan_inputs


@dataclass
class MaudInput:
    """Everything that is needed to run Maud."""

    config: MaudConfig
    kinetic_model: KineticModel
    experiments: List[Experiment]
    prior_input: PriorInput = PriorInput()
    init_input: InitInput = InitInput()
    parameters: ParameterSet = Field(init=False, exclude=True)
    stan_input_train: Dict = Field(init=False, exclude=True)
    stan_input_test: Dict = Field(init=False, exclude=True)

    def __post_init__(self):
        """Add attributes that depend on other ones."""
        self.parameters = get_maud_parameters(
            self.kinetic_model,
            self.experiments,
            self.prior_input,
            self.init_input,
        )
        self.stan_input_train, self.stan_input_test = get_stan_inputs(
            self.parameters,
            self.experiments,
            self.kinetic_model,
            self.config,
        )
        inits_dict = {}
        for p in map(
            lambda f: getattr(self.parameters, f.name),
            fields(self.parameters),
        ):
            inits_dict[p.name] = p.inits.inits_unscaled
            if p.inits.inits_scaled is not None:
                scaled_pref = "log_" if p.non_negative else ""
                inits_dict[scaled_pref + p.name + "_z"] = p.inits.inits_scaled
        self.inits_dict = inits_dict
