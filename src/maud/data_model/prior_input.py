"""Definitions of the user's input for priors. Directly read from toml."""

from typing import List, Optional, Union

from pydantic import root_validator, validator
from pydantic.dataclasses import dataclass


@dataclass
class IndPriorAtomInput:
    """Prior input for a single quantity."""

    metabolite: Optional[str] = None
    compartment: Optional[str] = None
    enzyme: Optional[str] = None
    reaction: Optional[str] = None
    experiment: Optional[str] = None
    phosphorylation_modifying_enzyme: Optional[str] = None
    modification_type: Optional[str] = None
    location: Optional[float] = None
    exploc: Optional[float] = None
    scale: Optional[float] = None
    pct1: Optional[float] = None
    pct99: Optional[float] = None

    @validator("scale")
    def scale_is_positive(cls, v):
        """Check that scale is positive."""

        if v is not None and v <= 0:
            raise ValueError("scale must be a positive number.")
        return v

    @root_validator(pre=False, skip_on_failure=True)
    def prior_is_specified_correctly(cls, values):
        """Check that location, scale etc are correct."""
        lc = values["location"]
        el = values["exploc"]
        sc = values["scale"]
        p1 = values["pct1"]
        p99 = values["pct99"]
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
        return values


@dataclass
class PriorMVNInput:
    """Prior input for a parameter with multivariate normal distribution."""

    ids: List[str]
    mean_vector: List[float]
    covariance_matrix: List[List[float]]


@dataclass
class PriorInput:
    """A full prior input."""

    dgf: Optional[Union[PriorMVNInput, List[IndPriorAtomInput]]] = None
    km: Optional[List[IndPriorAtomInput]] = None
    kcat: Optional[List[IndPriorAtomInput]] = None
    kcat_pme: Optional[List[IndPriorAtomInput]] = None
    ki: Optional[List[IndPriorAtomInput]] = None
    psi: Optional[List[IndPriorAtomInput]] = None
    dissociation_constant: Optional[List[IndPriorAtomInput]] = None
    transfer_constant: Optional[List[IndPriorAtomInput]] = None
    conc_unbalanced: Optional[List[IndPriorAtomInput]] = None
    drain: Optional[List[IndPriorAtomInput]] = None
    conc_enzyme: Optional[List[IndPriorAtomInput]] = None
    conc_pme: Optional[List[IndPriorAtomInput]] = None
