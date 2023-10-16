"""Provides model Experiment."""

from enum import Enum
from typing import List, Optional

from pydantic import BaseModel, Field, computed_field, field_validator

from maud.data_model.hardcoding import ID_SEPARATOR


class MeasurementType(str, Enum):
    """Possible types of measurement."""

    MIC = "mic"
    FLUX = "flux"
    ENZYME = "enzyme"


class Measurement(BaseModel):
    """Maud representation of a measurement."""

    experiment: str
    target_type: MeasurementType
    value: float = Field(kw_only=True, allow_inf_nan=False)
    error_scale: float = Field(kw_only=True, allow_inf_nan=False, gt=0)
    metabolite: Optional[str] = None
    compartment: Optional[str] = None
    reaction: Optional[str] = None
    enzyme: Optional[str] = None

    @computed_field
    def target_id(self) -> str:
        """Add target_id field."""
        if self.target_type == MeasurementType.MIC:
            assert self.metabolite is not None
            assert self.compartment is not None
            return ID_SEPARATOR.join([self.metabolite, self.compartment])
        elif self.target_type == MeasurementType.FLUX:
            assert self.reaction is not None
            return self.reaction
        else:
            assert self.enzyme is not None
            return self.enzyme


class EnzymeKnockout(BaseModel):
    """Maud representation of an enzyme being knocked out in an experiment."""

    experiment: str
    enzyme: str

    @computed_field
    def id(self) -> str:
        """Add id field."""
        return ID_SEPARATOR.join(["eko", self.experiment, self.enzyme])


class PhosphorylationModifyingEnzymeKnockout(BaseModel):
    """Maud representation of a pme being knocked out in an experiment."""

    experiment: str
    enzyme: str
    pme: str

    @computed_field
    def id(self) -> str:
        """Add id field."""
        return ID_SEPARATOR.join(["pko", self.experiment, self.enzyme])


class InitConcentration(BaseModel):
    """Indication of the initial value of a concentration in the ODE."""

    metabolite: str
    compartment: str
    value: float

    @computed_field
    def target_id(self) -> str:
        """Add target_id field."""
        return ID_SEPARATOR.join([self.metabolite, self.compartment])


class Experiment(BaseModel):
    """Maud representation of an experiment.

    This means a case where the boundary conditions and all measured quantities
    can be assumed to be the same - often the term "condition" is used for this.

    """

    id: str
    is_train: bool
    is_test: bool
    temperature: float = 298.15
    measurements: List[Measurement] = Field(default_factory=lambda: [])
    initial_state: List[InitConcentration] = Field(default_factory=lambda: [])
    enzyme_knockouts: List[EnzymeKnockout] = Field(default_factory=lambda: [])
    pme_knockouts: List[PhosphorylationModifyingEnzymeKnockout] = Field(
        default_factory=lambda: []
    )

    @field_validator("temperature")
    def temp_must_be_non_negative(cls, v):
        """Make sure the temperature isn't negative."""
        assert v >= 0
        return v


def parse_experiment(raw: dict):
    """Get an Experiment object from a dictionary that comes from toml.load."""
    special = {"measurements": [], "enzyme_knockouts": [], "pme_knockouts": []}
    not_special = {k: v for k, v in raw.items() if k not in special.keys()}
    for k, obj in zip(
        special.keys(),
        [Measurement, EnzymeKnockout, PhosphorylationModifyingEnzymeKnockout],
    ):
        if k in raw.keys():
            special[k] = [obj(experiment=raw["id"], **r) for r in raw[k]]
    return Experiment(**{**special, **not_special})
