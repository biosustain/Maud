"""Provides model MaudInput containing everything Maud needs to run."""

from typing import Dict, List

from pydantic import BaseModel, Field, computed_field

from maud.data_model.experiment import Experiment
from maud.data_model.kinetic_model import KineticModel
from maud.data_model.maud_config import MaudConfig
from maud.data_model.maud_init import InitInput
from maud.data_model.maud_parameter import ParameterSet
from maud.data_model.parameter_input import ParametersInput
from maud.getting_parameters import get_maud_parameters
from maud.getting_stan_inputs import get_stan_inputs


class MaudInput(BaseModel):
    """Everything that is needed to run Maud."""

    config: MaudConfig
    kinetic_model: KineticModel
    experiments: List[Experiment]
    parameters_input: ParametersInput = Field(default_factory=ParametersInput)
    init_input: InitInput = Field(default_factory=InitInput)

    @computed_field
    def parameters(self) -> ParameterSet:
        """Add the parameters field."""
        return get_maud_parameters(
            self.kinetic_model,
            self.experiments,
            self.parameters_input,
            self.init_input,
        )

    @computed_field
    def stan_input_train(self) -> Dict:
        """Add the stan_input_train field."""
        train, _ = get_stan_inputs(
            self.parameters,
            self.experiments,
            self.kinetic_model,
            self.config,
        )
        return train

    @computed_field
    def stan_input_test(self) -> Dict:
        """Add the stan_input_test field."""
        _, test = get_stan_inputs(
            self.parameters,
            self.experiments,
            self.kinetic_model,
            self.config,
        )
        return test

    @computed_field
    def inits_dict(self) -> Dict:
        """Add the inits_dict field."""
        inits_dict = {}
        for p in map(
            lambda k: getattr(self.parameters, k),
            self.parameters.__fields__.keys(),
        ):
            inits_dict[p.name] = p.inits.inits_unscaled
            if p.inits.inits_scaled is not None:
                scaled_pref = "log_" if p.non_negative else ""
                inits_dict[scaled_pref + p.name + "_z"] = p.inits.inits_scaled
        return inits_dict
