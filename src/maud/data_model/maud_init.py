"""Types that go in initial values specificaitons."""

from typing import List, Optional, Union

from pydantic.dataclasses import dataclass


@dataclass
class InitAtomInput:
    init: float
    metabolite: Optional[str] = None
    compartment: Optional[str] = None
    enzyme: Optional[str] = None
    reaction: Optional[str] = None
    experiment: Optional[str] = None


@dataclass
class Init1d:
    inits_unscaled: List[float]
    inits_scaled: Optional[List[float]] = None


@dataclass
class Init2d:
    inits_unscaled: List[List[float]]
    inits_scaled: Optional[List[List[float]]] = None


Init = Union[Init1d, Init2d]
ParamInitInput = Optional[List[InitAtomInput]]


@dataclass
class InitInput:
    dgf: ParamInitInput = None
    km: ParamInitInput = None
    kcat: ParamInitInput = None
    kcat_pme: ParamInitInput = None
    ki: ParamInitInput = None
    psi: ParamInitInput = None
    dissociation_constant: ParamInitInput = None
    transfer_constant: ParamInitInput = None
    conc_unbalanced: ParamInitInput = None
    drain: ParamInitInput = None
    conc_enzyme: ParamInitInput = None
    conc_pme: ParamInitInput = None
