"""Types that go in initial values specificaitons."""

from typing import List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from pydantic.dataclasses import dataclass

from maud.data_model.experiment import Measurement
from maud.data_model.hardcoding import ID_SEPARATOR
from maud.data_model.id_component import IdComponent
from maud.data_model.prior import IndPrior1d, IndPrior2d, PriorMVN


@dataclass
class InitAtomInput:
    """Maud representation of an init input for a single quantity."""

    init: float
    metabolite: Optional[str] = None
    compartment: Optional[str] = None
    enzyme: Optional[str] = None
    reaction: Optional[str] = None
    experiment: Optional[str] = None
    modification_type: Optional[str] = None


ParamInitInput = Optional[List[InitAtomInput]]


@dataclass
class InitInput:
    """Maud representation of a full init input."""

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


def get_init_atom_input_ids(
    iai: InitAtomInput,
    id_components: List[List[IdComponent]],
) -> Tuple[str, ...]:
    """Get an init's id from a userinput."""

    return tuple(
        ID_SEPARATOR.join([getattr(iai, c) for c in idci])
        for idci in id_components
    )


@dataclass(init=False)
class Init1d:
    """A 1 dimensional initial values specification."""

    inits_unscaled: List[float]
    inits_scaled: Optional[List[float]] = None

    def __init__(
        self,
        ids: List[List[str]],
        id_components: List[List[IdComponent]],
        prior: Union[IndPrior1d, PriorMVN],
        param_init_input: ParamInitInput,
        non_negative: bool,
        measurements: Optional[List[Measurement]] = None,
    ):
        inits_pd = pd.Series(prior.location, index=ids[0], dtype="float64")
        if non_negative:  # non-negative parameter location is on ln scale
            inits_pd = np.exp(inits_pd)
        if measurements is not None:
            for m in measurements:
                if m.target_id in inits_pd.index:
                    inits_pd.loc[m.target_id] = m.value
        if param_init_input is not None:
            for iai in param_init_input:
                iai_id = get_init_atom_input_ids(iai, id_components)[0]
                if iai_id in inits_pd.index:
                    inits_pd.loc[iai_id] = iai.init
        self.inits_unscaled = inits_pd.tolist()
        if isinstance(prior, PriorMVN):  # no need to rescale an MVN parameter
            self.inits_scaled = None
        else:
            inits_for_scaling = (
                np.log(inits_pd) if non_negative else inits_pd.copy()
            )
            loc_pd = pd.Series(prior.location, index=ids[0], dtype="float64")
            scale_pd = pd.Series(prior.scale, index=ids[0], dtype="float64")
            inits_pd_scaled = (inits_for_scaling - loc_pd) / scale_pd
            self.inits_scaled = inits_pd_scaled.tolist()


@dataclass(init=False)
class Init2d:
    """A 2 dimensional initial values specification."""

    inits_unscaled: List[List[float]]
    inits_scaled: Optional[List[List[float]]] = None

    def __init__(
        self,
        ids: List[List[str]],
        id_components: List[List[IdComponent]],
        prior: IndPrior2d,
        param_init_input: ParamInitInput,
        non_negative: bool,
        measurements: Optional[List[Measurement]] = None,
    ):
        inits_pd = pd.DataFrame(prior.location, index=ids[0], columns=ids[1])
        if non_negative:  # non-negative parameter location is on ln scale
            inits_pd = np.exp(inits_pd)
        if measurements is not None:
            for m in measurements:
                if (
                    m.target_id in inits_pd.columns
                    and m.experiment in inits_pd.index
                ):
                    inits_pd.loc[m.experiment, m.target_id] = m.value
        if param_init_input is not None:
            for iai in param_init_input:
                iai_id_row, iai_id_col = get_init_atom_input_ids(
                    iai, id_components
                )
                if (
                    iai_id_row in inits_pd.index
                    and iai_id_col in inits_pd.columns
                ):
                    inits_pd.loc[iai_id_row, iai_id_col] = iai.init
        inits_for_scaling = (
            np.log(inits_pd) if non_negative else inits_pd.copy()
        )
        loc_pd = pd.DataFrame(prior.location, index=ids[0], columns=ids[1])
        scale_pd = pd.DataFrame(prior.scale, index=ids[0], columns=ids[1])
        inits_pd_scaled = (inits_for_scaling - loc_pd) / scale_pd
        self.inits_unscaled = inits_pd.values.tolist()
        self.inits_scaled = inits_pd_scaled.values.tolist()


Init = Union[Init1d, Init2d]
