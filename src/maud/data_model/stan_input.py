"""Definitions of inputs for Maud's Stan models."""

from math import isnan
from typing import Dict, Sequence, Union

from pydantic import StrictFloat, StrictInt, validator
from pydantic.dataclasses import dataclass

from maud.utils import recursively_flatten_list


# Strict types are used here to stop pydantic from doing unwanted conversions.
# see these links for more discussion:
# - https://pydantic-docs.helpmanual.io/usage/models/#data-conversion
# - https://pydantic-docs.helpmanual.io/usage/types/#strict-types
StanDataValueBase = Union[StrictFloat, StrictInt]
StanDataValue = Union[
    StanDataValueBase,
    Sequence[StanDataValueBase],
    Sequence[Sequence[StanDataValueBase]],
    Sequence[Sequence[Sequence[StanDataValueBase]]],
]
StanInputDict = Dict[str, StanDataValue]


@dataclass
class StanInput:
    """Input for a Stan model."""

    stan_input_dict: StanInputDict

    @validator("stan_input_dict")
    def no_nans_allowed(cls, v):
        """Check that Stan input data isn't null and doesn't contain nulls."""
        for k, standata in v.items():
            msg = f"StanData {k} has non-numbers: {standata}"
            if isinstance(standata, Sequence):
                flat = recursively_flatten_list(standata)
                assert not any([isnan(x) for x in flat]), msg
            else:
                assert not isnan(standata), msg
        return v
