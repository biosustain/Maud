"""Functions for creating InitDict objects."""

from typing import List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from maud.data_model.hardcoding import ID_SEPARATOR
from maud.data_model.id_component import IdComponent
from maud.data_model.maud_init import Init1d, Init2d, InitAtomInput, ParamInitInput, Init
from maud.data_model.prior import Prior, PriorMVN
from maud.data_model.experiment import Measurement


def get_init_atom_input_ids(
    iai: InitAtomInput,
    id_components: List[List[IdComponent]],
) -> Tuple[str, ...]:
    """Get an init's id from a userinput."""
    return tuple(
        ID_SEPARATOR.join([getattr(iai, c) for c in idci])
        for idci in id_components
    )


def get_param_inits(
    ids: List[List[str]],
    id_components: List[List[IdComponent]],
    prior: Prior,
    init_input: ParamInitInput,
    non_negative: bool,
    measurements: Optional[List[Measurement]] = None
) -> Init:
    """Get initial values for a prarameter, possibly given a user input."""
    if len(ids) == 1:
        inits_pd = pd.Series(prior.location, index=ids[0])
        if non_negative:  # non-negative parameter location is on ln scale
            inits_pd = np.exp(inits_pd)
        if measurements is not None:
            for m in measurements:
                if m.target_id in inits_pd.index:
                    inits_pd.loc[m.target_id] = m.value
        if init_input is not None:
            for iai in init_input:
                init_id = get_init_atom_input_ids(iai, id_components)[0]
                inits_pd.loc[init_id] = iai.init
        if isinstance(prior, PriorMVN):  # no need to rescale an MVN parameter
            return Init1d(inits_pd.tolist())
        else:
            loc_trans = np.log(inits_pd) if non_negative else inits_pd.copy()
            scale_pd = pd.Series(prior.scale, index=ids[0])
            inits_pd_scaled = (loc_trans - loc_trans.mean()) / scale_pd
            return Init1d(inits_pd.tolist(), inits_pd_scaled.tolist())
    else:
        inits_pd = pd.DataFrame(prior.location, index=ids[0], columns=ids[1])
        if non_negative:  # non-negative parameter location is on ln scale
            inits_pd = np.exp(inits_pd)
        if measurements is not None:
            for m in measurements:
                if m.target_id in inits_pd.columns and m.experiment in inits_pd.index:
                    inits_pd.loc[m.experiment, m.target_id] = m.value
        if init_input is not None:
            for iai in init_input:
                init_id_row, init_id_col = get_init_atom_input_ids(
                    iai, id_components
                )
                inits_pd.loc[init_id_row, init_id_col] = iai.init
        loc_trans = np.log(inits_pd) if non_negative else inits_pd.copy()
        scale_pd = pd.DataFrame(prior.scale, index=ids[0], columns=ids[1])
        inits_pd_scaled = (loc_trans - loc_trans.mean()) / scale_pd
        return Init2d(
            inits_pd.values.tolist(), inits_pd_scaled.values.tolist()
        )


