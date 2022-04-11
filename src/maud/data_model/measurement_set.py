"""Provides dataclass MeasurementSet"""

from dataclasses import field
from enum import Enum
from typing import List, Optional

import pandas as pd
from pydantic.dataclasses import dataclass
from pydantic.main import BaseModel

from maud.data_model.hardcoding import ID_SEPARATOR


class MeasurementType(str, Enum):
    MIC = "mic"
    FLUX = "flux"
    ENZYME = "enzyme"


@dataclass
class Experiment:
    id: str
    is_train: bool
    is_test: bool


@dataclass
class EnzymeKnockout:
    experiment_id: str
    enzyme_id: str
    id: str = field(init=False)

    def __post_init__(self):
        self.id = ID_SEPARATOR.join(["eko", self.experiment_id, self.enzyme_id])


@dataclass
class PhosphorylationKnockout:
    experiment_id: str
    enzyme_id: str
    id: str = field(init=False)

    def __post_init__(self):
        self.id = ID_SEPARATOR.join(["pko", self.experiment_id, self.enzyme_id])


class MeasurementSet(BaseModel):
    """A complete set of measurements, including knockouts."""

    yconc: pd.DataFrame
    yflux: pd.DataFrame
    yenz: pd.DataFrame
    enzyme_knockouts: Optional[List[EnzymeKnockout]]
    phosphorylation_knockouts: Optional[List[PhosphorylationKnockout]]
    experiments: List[Experiment]

    class Config:
        arbitrary_types_allowed = True
