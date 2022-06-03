"""Provides dataclass MeasurementSet, representing a set of measurements."""

from enum import Enum
from typing import List, Optional

import pandas as pd
from pydantic import Field, validator
from pydantic.dataclasses import dataclass

from maud.data_model.hardcoding import ID_SEPARATOR


class MSConfig:
    """Config for MeasurementSet, allowing it to contain pandas objects."""

    arbitrary_types_allowed = True


class MeasurementType(str, Enum):
    """Possible types of measurement."""

    MIC = "mic"
    FLUX = "flux"
    ENZYME = "enzyme"


@dataclass
class Experiment:
    """Maud representation of an experiment.

    This means a case where the boundary conditions and all measured quantities
    can be assumed to be the same - often the term "condition" is used for this.

    """

    id: str
    is_train: bool
    is_test: bool
    temperature: float = 298.15

    def __post_init__(self):
        """Add set temperature to default if it is None."""
        if self.temperature is None:
            self.temperature = 298.15

    @validator("temperature")
    def temp_must_be_non_negative(cls, v):
        """Make sure the temperature isn't negative."""
        assert v >= 0
        return v


@dataclass
class EnzymeKnockout:
    """Maud representation of an enzyme being knocked out in an experiment."""

    experiment_id: str
    enzyme_id: str
    id: str = Field(init=False, exclude=True)

    def __post_init__(self):
        """Add id field."""
        self.id = ID_SEPARATOR.join(
            ["eko", self.experiment_id, self.enzyme_id]
        )


@dataclass
class PhosphorylationModifyingEnzymeKnockout:
    """Maud representation of a pme being knocked out in an experiment."""

    experiment_id: str
    pme_id: str
    id: str = Field(init=False, exclude=True)

    def __post_init__(self):
        """Add id field."""
        self.id = ID_SEPARATOR.join(
            ["pko", self.experiment_id, self.enzyme_id]
        )


@dataclass(config=MSConfig)
class MeasurementSet:
    """A complete set of measurements, including knockouts."""

    yconc: pd.DataFrame
    yflux: pd.DataFrame
    yenz: pd.DataFrame
    enzyme_knockouts: Optional[List[EnzymeKnockout]]
    pme_knockouts: Optional[List[PhosphorylationModifyingEnzymeKnockout]]
    experiments: List[Experiment]
