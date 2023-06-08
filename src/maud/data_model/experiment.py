"""Provides dataclass Experiment."""

from enum import Enum
from typing import List, Optional

from pydantic import Field, root_validator, validator
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
class Measurement:
    """Maud representation of a measurement."""

    experiment: str
    target_type: MeasurementType
    value: float = Field(kw_only=True, allow_inf_nan=False)
    error_scale: Optional[float] = Field(
        kw_only=True, allow_inf_nan=False, gt=0, default=None
    )
    metabolite: Optional[str] = None
    compartment: Optional[str] = None
    reaction: Optional[str] = None
    enzyme: Optional[str] = None
    target_id: str = Field(default=None, init=False, exclude=True)
    # if False, this is used only for specifying an initial value
    likelihood: bool = Field(default=True)

    def __post_init__(self):
        """Add target_id field."""
        if self.target_type == MeasurementType.MIC:
            self.target_id = ID_SEPARATOR.join(
                [self.metabolite, self.compartment]
            )
        elif self.target_type == MeasurementType.FLUX:
            self.target_id = self.reaction
        elif self.target_type == MeasurementType.ENZYME:
            self.target_id = self.enzyme

    @root_validator(pre=False, skip_on_failure=True)
    def error_scale_must_exist_if_likelihood(cls, values):
        """Verify the right combinations of error_scale and likelihood."""
        assert values["reaction"] is None or values["likelihood"], (
            "likelihood=False is not implemented for reaction types"
            f"({values['reaction']}, {values['experiment']})"
        )
        ident = (
            values["reaction"]
            if values["reaction"] is not None
            else values["metabolite"]
        )
        assert values["error_scale"] is not None or not values["likelihood"], (
            f"{ident} measurement (exp={values['experiment']}) must provide"
            " error_scale if used for likelihood"
        )
        assert not (
            values["error_scale"] is not None and not values["likelihood"]
        ), (
            f"{ident} measurement (exp={values['experiment']}) has an"
            " error_scale but likelihood is false"
        )
        return values


@dataclass
class EnzymeKnockout:
    """Maud representation of an enzyme being knocked out in an experiment."""

    experiment: str
    enzyme: str
    id: str = Field(init=False, exclude=True)

    def __post_init__(self):
        """Add id field."""
        self.id = ID_SEPARATOR.join(["eko", self.experiment, self.enzyme])


@dataclass
class PhosphorylationModifyingEnzymeKnockout:
    """Maud representation of a pme being knocked out in an experiment."""

    experiment: str
    pme: str
    id: str = Field(init=False, exclude=True)

    def __post_init__(self):
        """Add id field."""
        self.id = ID_SEPARATOR.join(["pko", self.experiment, self.enzyme])


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
    measurements: List[Measurement] = Field(default_factory=lambda: [])
    enzyme_knockouts: List[EnzymeKnockout] = Field(default_factory=lambda: [])
    pme_knockouts: List[PhosphorylationModifyingEnzymeKnockout] = Field(
        default_factory=lambda: []
    )

    @validator("temperature")
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
