"""Provides model MaudInput containing everything Maud needs to run."""

from typing import Dict, List

from pydantic import BaseModel, Field, computed_field

from maud.data_model.experiment import Experiment
from maud.data_model.kinetic_model import KineticModel
from maud.data_model.maud_config import MaudConfig
from maud.data_model.maud_init import InitInput
from maud.data_model.maud_parameter import MaudParameter
from maud.data_model.parameter_input import ParameterSetInput
from maud.data_model.parameter_set import ParameterSet
from maud.getting_stan_inputs import get_stan_inputs


class MaudInput(BaseModel):
    """Everything that is needed to run Maud."""

    config: MaudConfig
    kinetic_model: KineticModel
    experiments: List[Experiment]
    parameter_set_input: ParameterSetInput = Field(
        default_factory=ParameterSetInput
    )
    init_input: InitInput = Field(default_factory=InitInput)

    @computed_field
    def parameters(self) -> ParameterSet:
        """Add the parameters field."""
        return ParameterSet(
            kinetic_model=self.kinetic_model,
            experiments=self.experiments,
            parameter_set_input=self.parameter_set_input,
            init_input=self.init_input,
        )

    @computed_field
    def stan_input_train(self) -> Dict:
        """Add the stan_input_train field."""
        train, _ = get_stan_inputs(
            parameters=self.parameters,
            experiments=self.experiments,
            kinetic_model=self.kinetic_model,
            config=self.config,
        )
        return train

    @computed_field
    def stan_input_test(self) -> Dict:
        """Add the stan_input_test field."""
        _, test = get_stan_inputs(
            parameters=self.parameters,
            experiments=self.experiments,
            kinetic_model=self.kinetic_model,
            config=self.config,
        )
        return test

    @computed_field
    def inits_dict(self) -> Dict:
        """Add the inits_dict field."""
        inits_dict = {}
        params = [
            getattr(self.parameters, p)
            for p in self.parameters.dict().keys()
            if isinstance(getattr(self.parameters, p), MaudParameter)
        ]
        for param in params:
            inits_dict[param.name] = param.inits.inits_unscaled
            if param.inits.inits_scaled is not None:
                scaled_pref = "log_" if param.non_negative else ""
                inits_dict[
                    scaled_pref + param.name + "_z"
                ] = param.inits.inits_scaled
            if param.fixed_ids is not None:
                met_to_init = dict(
                    zip(param.inits.ids[0], param.inits.inits_unscaled)
                )
                inits_dict[param.name + "_free"] = [
                    init
                    for met, init in met_to_init.items()
                    if met not in param.fixed_ids[0]
                ]
        return inits_dict
