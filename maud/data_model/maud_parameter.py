"""Provides model MaudParameter."""

from typing import List, Optional, Union

from pydantic import BaseModel, computed_field, field_validator, model_validator

from maud.data_model.experiment import Measurement
from maud.data_model.id_component import IdComponent
from maud.data_model.maud_init import Init, Init1d, Init2d, InitAtomInput
from maud.data_model.parameter_input import (
    ParameterInputAtom,
    ParameterInputMVN,
)
from maud.data_model.prior import IndPrior1d, IndPrior2d, PriorMVN


class MaudParameter(BaseModel):
    """A parameter in Maud's statistical model."""

    name: str
    shape_names: List[str]
    id_components: List[List[IdComponent]]
    non_negative: bool
    default_scale: float
    default_loc: float
    prior_in_test_model: bool
    prior_in_train_model: bool
    user_input: Optional[Union[List[ParameterInputAtom], ParameterInputMVN]]
    init_input: Optional[List[InitAtomInput]]
    ids: List[List[str]]
    split_ids: List[List[List[str]]]
    measurements: Optional[List[Measurement]] = None

    @computed_field
    def prior(self) -> Union[IndPrior1d, IndPrior2d, PriorMVN]:
        """Return a prior, calculated from the user input."""
        if self.name == "dgf":
            initialiser = PriorMVN
        elif len(self.shape_names) == 1:
            initialiser = IndPrior1d
        else:
            initialiser = IndPrior2d
        return initialiser(
            user_input=self.user_input,
            ids=self.ids,
            id_components=self.id_components,
            non_negative=self.non_negative,
            default_loc=self.default_loc,
            default_scale=self.default_scale,
        )

    @computed_field
    def inits(self) -> Init:
        """Add the inits field."""
        initialiser = Init1d if len(self.shape_names) == 1 else Init2d
        return initialiser(
            ids=self.ids,
            id_components=self.id_components,
            prior=self.prior,
            param_init_input=self.init_input,
            non_negative=self.non_negative,
            measurements=self.measurements,
        )

    @field_validator("id_components")
    def id_components_have_same_length(cls, v):
        """Make sure that id_components contains lists with the same length."""
        first_length = len(v[0])
        assert all([len(x) == first_length for x in v])
        return v

    @model_validator(mode="after")
    def split_ids_exist_if_needed(self):
        """Check split ids exist when there are non-trivial id components."""
        if any(len(idc) > 1 for idc in self.id_components):
            assert self.split_ids is not None
        return self
