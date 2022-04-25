"""Provides dataclass MeasurementSet, representing a set of measurements."""

from dataclasses import field
from enum import Enum
from typing import List, Optional

import pandas as pd
from pydantic.dataclasses import dataclass
from pydantic.main import BaseModel

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


@dataclass
class EnzymeKnockout:
    """Maud representation of an enzyme being knocked out in an experiment."""

    experiment_id: str
    enzyme_id: str
    id: str = field(init=False)

    def __post_init__(self):
        """Add id field."""
        self.id = ID_SEPARATOR.join(
            ["eko", self.experiment_id, self.enzyme_id]
        )


@dataclass
class PhosphorylationKnockout:
    """Maud representation of a phosphorylation being prevented in an experiment."""

    experiment_id: str
    enzyme_id: str
    id: str = field(init=False)

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
    phosphorylation_knockouts: Optional[List[PhosphorylationKnockout]]
    experiments: List[Experiment]
