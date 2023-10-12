"""Definitions of the user's input for priors. Directly read from toml."""

from typing import Dict, List, Optional, Union

from pydantic import (
    BaseModel,
    Field,
    PositiveFloat,
    field_validator,
    model_validator,
)


class ParameterInputAtom(BaseModel):
    """Parameter input for a single quantity."""

    metabolite: Optional[str] = None
    compartment: Optional[str] = None
    enzyme: Optional[str] = None
    reaction: Optional[str] = None
    experiment: Optional[str] = None
    phosphorylation_modifying_enzyme: Optional[str] = None
    modification_type: Optional[str] = None
    location: Optional[float] = None
    exploc: Optional[float] = None
    scale: Optional[PositiveFloat] = None
    pct1: Optional[float] = None
    pct99: Optional[float] = None
    init: Optional[float] = None
    fixed_value: Optional[float] = None

    @field_validator("scale")
    def scale_is_positive(cls, v):
        """Check that scale is positive."""
        if v is not None and v <= 0:
            raise ValueError("scale must be a positive number.")
        return v

    @model_validator(mode="before")
    def prior_is_specified_correctly(cls, data):
        """Check that location, scale etc are correct."""
        lc = data["location"] if "location" in data.keys() else None
        el = data["exploc"] if "exploc" in data.keys() else None
        sc = data["scale"] if "scale" in data.keys() else None
        p1 = data["pct1"] if "pct1" in data.keys() else None
        p99 = data["pct99"] if "pct99" in data.keys() else None
        happy_cases = [
            {"not_none": [lc, sc], "none": [el, p1, p99]},
            {"not_none": [el, sc], "none": [lc, p1, p99]},
            {"not_none": [p1, p99], "none": [lc, el, sc]},
        ]
        good = [
            all(v is not None for v in case["not_none"])
            and all(v is None for v in case["none"])
            for case in happy_cases
        ]
        if not any(good):
            raise ValueError(
                "Set one out of the following pairs of attributes: "
                "location and scale, exploc and scale, or pct1 and pct99."
            )
        return data


class ParameterInputMVN(BaseModel):
    """User input for a parameter with multivariate normal prior."""

    ids: List[str]
    fixed_values: Optional[Dict[str, float]] = None
    mean_vector: List[float]
    covariance_matrix: List[List[float]]


class ParameterSetInput(BaseModel):
    """User input for all parameters."""

    dgf: Optional[Union[ParameterInputMVN, List[ParameterInputAtom]]] = Field(
        default_factory=list
    )
    km: Optional[List[ParameterInputAtom]] = Field(default_factory=list)
    kcat: Optional[List[ParameterInputAtom]] = Field(default_factory=list)
    kcat_pme: Optional[List[ParameterInputAtom]] = Field(default_factory=list)
    ki: Optional[List[ParameterInputAtom]] = Field(default_factory=list)
    psi: Optional[List[ParameterInputAtom]] = Field(default_factory=list)
    dissociation_constant: Optional[List[ParameterInputAtom]] = Field(
        default_factory=list
    )
    transfer_constant: Optional[List[ParameterInputAtom]] = Field(
        default_factory=list
    )
    conc_unbalanced: Optional[List[ParameterInputAtom]] = Field(
        default_factory=list
    )
    drain: Optional[List[ParameterInputAtom]] = Field(default_factory=list)
    conc_enzyme: Optional[List[ParameterInputAtom]] = Field(
        default_factory=list
    )
    conc_pme: Optional[List[ParameterInputAtom]] = Field(default_factory=list)
