"""Types that go in initial values specificaitons."""

from typing import List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from pydantic import BaseModel, computed_field

from maud.data_model.experiment import Measurement
from maud.data_model.hardcoding import ID_SEPARATOR
from maud.data_model.id_component import IdComponent
from maud.data_model.prior import IndPrior1d, IndPrior2d, PriorMVN


class InitAtomInput(BaseModel):
    """Maud representation of an init input for a single quantity."""

    init: float
    metabolite: Optional[str] = None
    compartment: Optional[str] = None
    enzyme: Optional[str] = None
    reaction: Optional[str] = None
    experiment: Optional[str] = None
    modification_type: Optional[str] = None


ParamInitInput = Optional[List[InitAtomInput]]


class InitInput(BaseModel):
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


class Init1d(BaseModel):
    """A 1 dimensional initial values specification."""

    ids: List[List[str]]
    id_components: List[List[IdComponent]]
    prior: Union[IndPrior1d, PriorMVN]
    param_init_input: ParamInitInput
    non_negative: bool
    measurements: Optional[List[Measurement]] = None

    @computed_field
    def inits_unscaled(self) -> List[float]:
        """Add the inits_unscaled field."""
        inits_pd = pd.Series(
            self.prior.location, index=self.ids[0], dtype="float64"
        )
        if self.non_negative:  # non-negative parameter location is on ln scale
            inits_pd = np.exp(inits_pd)
        if self.measurements is not None:
            for m in self.measurements:
                if m.target_id in inits_pd.index:
                    inits_pd.loc[m.target_id] = m.value
        if self.param_init_input is not None:
            for iai in self.param_init_input:
                iai_id = get_init_atom_input_ids(iai, self.id_components)[0]
                if iai_id in inits_pd.index:
                    inits_pd.loc[iai_id] = iai.init
        return inits_pd.tolist()

    @computed_field
    def inits_scaled(self) -> Optional[List[float]]:
        """Add the inits_scaled field."""
        if isinstance(
            self.prior, PriorMVN
        ):  # no need to rescale an MVN parameter
            return None
        inits_pd = pd.Series(
            self.inits_unscaled, index=self.ids[0], dtype="float64"
        )
        inits_for_scaling = (
            np.log(inits_pd) if self.non_negative else inits_pd.copy()
        )
        loc_pd = pd.Series(
            self.prior.location, index=self.ids[0], dtype="float64"
        )
        scale_pd = pd.Series(
            self.prior.scale, index=self.ids[0], dtype="float64"
        )
        inits_pd_scaled = (inits_for_scaling - loc_pd) / scale_pd
        return inits_pd_scaled.tolist()


class Init2d(BaseModel):
    """A 2 dimensional initial values specification."""

    ids: List[List[str]]
    id_components: List[List[IdComponent]]
    prior: IndPrior2d
    param_init_input: ParamInitInput
    non_negative: bool
    measurements: Optional[List[Measurement]] = None

    @computed_field
    def inits_unscaled(self) -> List[List[float]]:
        """Add the inits_unscaled field."""
        inits_pd = pd.DataFrame(
            self.prior.location, index=self.ids[0], columns=self.ids[1]
        )
        if self.non_negative:  # non-negative parameter location is on ln scale
            inits_pd = np.exp(inits_pd)
        if self.measurements is not None:
            for m in self.measurements:
                if (
                    m.target_id in inits_pd.columns
                    and m.experiment in inits_pd.index
                ):
                    inits_pd.loc[m.experiment, m.target_id] = m.value
        if self.param_init_input is not None:
            for iai in self.param_init_input:
                iai_id_row, iai_id_col = get_init_atom_input_ids(
                    iai, self.id_components
                )
                if (
                    iai_id_row in inits_pd.index
                    and iai_id_col in inits_pd.columns
                ):
                    inits_pd.loc[iai_id_row, iai_id_col] = iai.init
        return inits_pd.values.tolist()

    @computed_field
    def inits_scaled(self) -> Optional[List[List[float]]]:
        """Add the inits_scaled field."""
        inits_pd = pd.DataFrame(
            self.inits_unscaled, index=self.ids[0], columns=self.ids[1]
        )
        inits_for_scaling = (
            np.log(inits_pd) if self.non_negative else inits_pd.copy()
        )
        loc_pd = pd.DataFrame(
            self.prior.location, index=self.ids[0], columns=self.ids[1]
        )
        scale_pd = pd.DataFrame(
            self.prior.scale, index=self.ids[0], columns=self.ids[1]
        )
        inits_pd_scaled = (inits_for_scaling - loc_pd) / scale_pd
        return inits_pd_scaled.values.tolist()


Init = Union[Init1d, Init2d]
