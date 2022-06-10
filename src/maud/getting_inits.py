"""Functions for creating InitDict objects."""

from typing import Optional, Union

import numpy as np
import pandas as pd

from maud.data_model.hardcoding import ID_SEPARATOR
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
    inits = {}
    for k, p in priors.__dict__.items():
        if k.startswith("__"):
            continue
        inits[p.stan_variable.name] = get_inits_for_param(user_inits, p)
    return rescale_inits(inits, priors)


def get_inits_for_param(
    user_init_df: Optional[pd.DataFrame],
    prior: Union[IndPrior1d, IndPrior2d, MultiVariateNormalPrior1d],
) -> Union[pd.Series, pd.DataFrame]:
    """Process a parameter and a dataframe of user-provided initial values."""
    inits_pd = (
        np.exp(prior.location).copy()
        if prior.stan_variable.non_negative
        else prior.location
    )
    if user_init_df is None:
        assert isinstance(inits_pd, pd.Series) or isinstance(
            inits_pd, pd.DataFrame
        )
        return inits_pd
    user = user_init_df.loc[
        lambda df: df["parameter"] == prior.stan_variable.name
    ]
    if len(prior.location) == 0:
        return pd.Series(prior.location).astype(float)
    elif isinstance(inits_pd, pd.Series):
        id_components = [c.value for c in prior.stan_variable.id_components[0]]
        if len(user) > 0:
            for i, row in user.iterrows():
                param_id = ID_SEPARATOR.join([row[c] for c in id_components])
                inits_pd.loc[param_id] = row["value"]
        return inits_pd
    elif isinstance(inits_pd, pd.DataFrame):
        id_components = [c.value for c in prior.stan_variable.id_components[1]]
        if len(user) > 0:
            for i, row in user.iterrows():
                non_experiment_id = ID_SEPARATOR.join(
                    [row[c] for c in id_components]
                )
                inits_pd.loc[row["experiment"], non_experiment_id] = row[
                    "value"
                ]
        return inits_pd
    else:
        raise ValueError("Encountered non-pandas inits: " + str(inits_pd))


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
