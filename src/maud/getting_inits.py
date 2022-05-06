"""Functions for creating InitDict objects."""

from typing import Optional, Union

import numpy as np
import pandas as pd

from maud.data_model.maud_init import InitDict
from maud.data_model.prior_set import (
    IndPrior1d,
    IndPrior2d,
    MultiVariateNormalPrior1d,
    PriorSet,
)


def get_inits(
    priors: PriorSet, user_inits: Optional[pd.DataFrame]
) -> InitDict:
    """Get a dictionary of initial values.

    :param priors: Priorset object

    :param user_inits_path: path to a csv of user-specified initial parameter
    values

    """

    prior_inits = {
        p.stan_variable.name: np.exp(p.location)
        if p.stan_variable.non_negative
        else p.location
        for k, p in priors.__dict__.items()
        if not k.startswith("__")
    }
    inits = prior_inits.copy()
    if user_inits is not None:
        for k, p in priors.__dict__.items():
            if k.startswith("__"):
                continue
            if p.location.empty:
                continue
            user = get_user_inits_for_param(user_inits, p)
            default = prior_inits[p.stan_variable.name]
            if isinstance(prior_inits[p.stan_variable.name], pd.DataFrame):
                default = default.stack()
            combined = pd.Series(
                np.where(
                    user.reindex(default.index).notnull(),
                    user.reindex(default.index),
                    default,
                ),
                index=default.index,
            )
            if isinstance(prior_inits[p.stan_variable.name], pd.DataFrame):
                inits[p.stan_variable.name] = combined.unstack()
            else:
                inits[p.stan_variable.name] = combined
    return rescale_inits(inits, priors)


def get_user_inits_for_param(
    u: pd.DataFrame,
    p: Union[IndPrior1d, IndPrior2d, MultiVariateNormalPrior1d],
) -> pd.Series:
    """Process a parameter and a dataframe of user-provided initial values."""
    if len(p.location) == 0:
        v = pd.Series(p.location).astype(float).values
        assert isinstance(v, np.ndarray)
        return v.tolist()
    elif isinstance(p, IndPrior1d) or isinstance(p, MultiVariateNormalPrior1d):
        v = (
            u.loc[lambda df: df["parameter_name"] == p.stan_variable.name]
            .set_index(p.location.index.names)["value"]
            .astype(float)
            .values
        )
        assert isinstance(v, np.ndarray)
        return v.tolist()
    elif isinstance(p, IndPrior2d):
        v = (
            u.loc[lambda df: df["parameter_name"] == p.stan_variable.name]
            .set_index([p.location.index.name, p.location.columns.name])[
                "value"
            ]
            .astype(float)
            .values
        )
        assert isinstance(v, np.ndarray)
        return v.tolist()
    else:
        raise ValueError("Unrecognised prior type: " + str(type(p)))


def rescale_inits(inits: dict, priors: PriorSet) -> InitDict:
    """Augment a dictionary of inits with equivalent normalised values.

    :param inits: original inits

    :param priors: PriorSet object used to do the normalising
    """
    rescaled = {}
    for (n, i), prior in zip(inits.items(), priors.__dict__.values()):
        rescaled_name = (
            f"log_{n}_z" if prior.stan_variable.non_negative else f"{n}_z"
        )
        if isinstance(prior, MultiVariateNormalPrior1d):
            continue
        else:
            loc_current = np.log(i) if prior.stan_variable.non_negative else i
            rescaled[rescaled_name] = (
                ((loc_current - prior.location) / prior.scale)
                .astype(float)
                .values.tolist()
            )
    return {
        **{k: v.astype(float).values.tolist() for k, v in inits.items()},
        **rescaled,
    }
