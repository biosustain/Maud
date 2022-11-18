"""Definitions of the user's input for priors. Directly read from toml."""

from pydantic.dataclasses import dataclass
from pydantic import validator, root_validator
from typing import Optional, List


@dataclass
class PriorInput0d:
    metabolite: Optional[str] = None
    compartment: Optional[str] = None
    enzyme: Optional[str] = None
    reaction: Optional[str] = None
    experiment: Optional[str] = None
    exploc: Optional[float] = None
    scale: Optional[float] = None
    pct1: Optional[float] = None
    pct99: Optional[float] = None

    @validator("scale")
    def scale_is_positive(cls, v):
        if v <= 0:
            raise ValueError("scale must be a positive number.")
        return v

    @root_validator
    def prior_is_specified_correctly(cls, values):
        el = values["exploc"]
        sc = values["scale"]
        p1 = values["pct1"]
        p99 = values["pct99"]
        good = [
            all([el is not None, sc is not None, p1 is None, p99 is None]),
            all([el is None, sc is None, p1 is not None, p99 is not None]),
        ]
        if not any(good):
            raise ValueError("Either set exploc and scale, or pct1 and pct99.")
        return values


@dataclass
class PriorInputMVN:
    ids: List[str]
    mean_vector: List[float]
    covariance_matrix: List[List[float]]


@dataclass
class PriorInput:
    dgf: Optional[PriorInputMVN] = None
    km: Optional[List[PriorInput0d]] = None
    kcat: Optional[List[PriorInput0d]] = None
    kcat_pme: Optional[List[PriorInput0d]] = None
    ki: Optional[List[PriorInput0d]] = None
    psi: Optional[List[PriorInput0d]] = None
    dissociation_constant: Optional[List[PriorInput0d]] = None
    transfer_constant: Optional[List[PriorInput0d]] = None
    conc_unbalanced: Optional[List[PriorInput0d]] = None
    drain: Optional[List[PriorInput0d]] = None
    conc_enzyme: Optional[List[PriorInput0d]] = None
    conc_pme: Optional[List[PriorInput0d]] = None
