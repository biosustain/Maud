"""Functions for creating prior objects from PriorInput objects.

This module handles setting priors to default values and assigning priors
consistently with their ids.

"""

from typing import List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from maud.data_model.hardcoding import ID_SEPARATOR
from maud.data_model.maud_parameter import IdComponent
from maud.data_model.prior import IndPrior1d, IndPrior2d, PriorMVN
from maud.data_model.prior_input import IndPriorAtomInput, PriorMVNInput
from maud.utils import (
    get_lognormal_parameters_from_quantiles,
    get_normal_parameters_from_quantiles,
)


def get_loc_and_scale(
    ipai: IndPriorAtomInput, non_negative: bool
) -> Tuple[float, float]:
    """Get location and scale from a user input prior.

    This is non-trivial because of the "exploc" input and the option to input
    priors as percentiles.

    """
    if ipai.location is not None and ipai.scale is not None:
        return ipai.location, ipai.scale
    elif ipai.exploc is not None and ipai.scale is not None:
        return np.log(ipai.exploc), ipai.scale
    elif ipai.pct1 is not None and ipai.pct99 is not None:
        quantfunc = (
            get_lognormal_parameters_from_quantiles
            if non_negative
            else get_normal_parameters_from_quantiles
        )
        return quantfunc(ipai.pct1, 0.01, ipai.pct99, 0.99)
    else:
        raise ValueError("Incorrect prior input.")


def unpack_ind_prior_atom_input(
    ipai: IndPriorAtomInput,
    id_components: List[List[IdComponent]],
    non_negative: bool,
) -> Tuple[List[str], float, float]:
    """Get a 0d prior from a user input."""
    ids = [
        ID_SEPARATOR.join([getattr(ipai, c) for c in idci])
        for idci in id_components
    ]
    loc, scale = get_loc_and_scale(ipai, non_negative)
    return ids, loc, scale


def get_ind_prior_1d(
    pi: Optional[List[IndPriorAtomInput]],
    ids: List[str],
    id_components: List[List[IdComponent]],
    non_negative: bool,
    default_loc: float,
    default_scale: float,
) -> IndPrior1d:
    """Get an independent 1d prior from a prior input and StanVariable."""
    if len(ids) == 0:
        return IndPrior1d(location=[], scale=[])
    loc_series = pd.Series(default_loc, index=ids)
    scale_series = pd.Series(default_scale, index=ids)
    if pi is not None:
        for ipai in pi:
            ids_i, loc_i, scale_i = unpack_ind_prior_atom_input(
                ipai, id_components, non_negative
            )
            if ids_i[0] in loc_series.index:
                loc_series.update({ids_i[0]: loc_i})
                scale_series.update({ids_i[0]: scale_i})
    return IndPrior1d(loc_series.tolist(), scale_series.tolist())


def get_ind_prior_2d(
    pi: Optional[List[IndPriorAtomInput]],
    ids: List[List[str]],
    id_components: List[List[IdComponent]],
    non_negative: bool,
    default_loc: float,
    default_scale: float,
) -> IndPrior2d:
    """Get an independent 2d prior from a prior input and StanVariable."""
    if any(len(ids_i) == 0 for ids_i in ids):
        return IndPrior2d(location=[[]], scale=[[]])
    loc_df = pd.DataFrame(float(default_loc), index=ids[0], columns=ids[1])
    scale_df = pd.DataFrame(float(default_scale), index=ids[0], columns=ids[1])
    if pi is not None:
        for ipai in pi:
            ids_i, loc_i, scale_i = unpack_ind_prior_atom_input(
                ipai, id_components, non_negative
            )
            if ids_i[0] in loc_df.index and ids_i[1] in loc_df.columns:
                loc_df.loc[ids_i[0], ids_i[1]] = loc_i
                scale_df.loc[ids_i[0], ids_i[1]] = scale_i
    return IndPrior2d(loc_df.values.tolist(), scale_df.values.tolist())


def get_mvn_prior(
    pi: Optional[Union[List[IndPriorAtomInput], PriorMVNInput]],
    ids: List[str],
    id_components: List[List[IdComponent]],
    non_negative: bool,
    default_loc: float,
    default_scale: float,
) -> PriorMVN:
    """Get a multivariate normal prior from a prior input and StanVariable."""
    loc_series = pd.Series(default_loc, index=ids)
    cov_df = pd.DataFrame(
        np.diagflat(np.tile(default_scale, len(ids))), index=ids, columns=ids
    )
    if isinstance(pi, PriorMVNInput):
        loc_series = pd.Series(pi.mean_vector, index=pi.ids).reindex(ids)
        cov_df = (
            pd.DataFrame(pi.covariance_matrix, index=pi.ids, columns=pi.ids)
            .reindex(ids)
            .reindex(columns=ids)
        )
    elif isinstance(pi, list):
        for ipai in pi:
            ids_i, loc_i, cov_ii = unpack_ind_prior_atom_input(
                ipai, id_components, non_negative
            )
            loc_series.loc[ids_i[0]] = loc_i
            cov_df.loc[ids_i[0], ids_i[0]] = cov_ii
    return PriorMVN(loc_series.tolist(), cov_df.values.tolist())
