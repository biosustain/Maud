"""Functions for creating *Prior* and PriorSet objects."""

from functools import partial

import numpy as np
import pandas as pd

from maud.data_model.hardcoding import ID_SEPARATOR
from maud.data_model.prior_set import (
    IndPrior1d,
    IndPrior2d,
    MultiVariateNormalPrior1d,
    PriorSet,
    UserPriorInput,
)
from maud.data_model.stan_variable_set import StanVariable, StanVariableSet
from maud.utils import (
    get_lognormal_parameters_from_quantiles,
    get_normal_parameters_from_quantiles,
    series_to_diag_df,
)


def load_1d_prior(
    user_df: pd.DataFrame, stan_variable: StanVariable
) -> IndPrior1d:
    """Get an IndPrior1d object from a dataframe and StanVariable.

    The StanVariable provides defaults and is included in the return value.
    """
    location = pd.Series(stan_variable.default_loc, index=stan_variable.ids[0])
    scale = pd.Series(stan_variable.default_scale, index=stan_variable.ids[0])
    user_df_sv = user_df.loc[lambda df: df["parameter"] == stan_variable.name]
    id_cols = [idc.value for idc in stan_variable.id_components[0]]
    if all(c not in user_df_sv.columns for c in id_cols):
        return IndPrior1d(stan_variable, location, scale)
    ids = user_df_sv[id_cols].apply(ID_SEPARATOR.join, axis=1).tolist()
    qf = (
        partial(get_lognormal_parameters_from_quantiles, p1=0.01, p2=0.99)
        if stan_variable.non_negative
        else partial(get_normal_parameters_from_quantiles, p1=0.01, p2=0.99)
    )
    pct_loc, pct_scale = qf(x1=user_df_sv["pct1"], x2=user_df_sv["pct99"])
    ls_loc, ls_scale = user_df_sv["location"], user_df_sv["scale"]
    # this is because prior locations for non-negative variables are enterred
    # unlogged:
    if stan_variable.non_negative:
        ls_loc = np.log(ls_loc)
    is_ls = user_df_sv[["location", "scale"]].notnull().all(axis=1)
    location.loc[ids] = np.where(is_ls, ls_loc, pct_loc)
    scale.loc[ids] = np.where(is_ls, ls_scale, pct_scale)
    return IndPrior1d(stan_variable, location, scale)


def load_2d_prior(
    user_df: pd.DataFrame, stan_variable: StanVariable
) -> IndPrior2d:
    """Get an IndPrior2d object from a dataframe and StanVariable.

    The StanVariable provides defaults and is included in the return value.
    """
    loc = pd.DataFrame(
        stan_variable.default_loc,
        index=stan_variable.ids[0],
        columns=stan_variable.ids[1],
    )
    scale = pd.DataFrame(
        stan_variable.default_scale,
        index=stan_variable.ids[0],
        columns=stan_variable.ids[1],
    )
    user = user_df.loc[lambda df: df["parameter"] == stan_variable.name].copy()
    row_id_cols = [idc.value for idc in stan_variable.id_components[1]]
    if all(c not in user.columns for c in row_id_cols):
        return IndPrior2d(stan_variable, loc, scale)
    user["col_id"] = (
        pd.NA
        if user[row_id_cols].empty
        else (user[row_id_cols].apply(ID_SEPARATOR.join, axis=1))
    )
    user["row_id"] = pd.NA if user["experiment"].empty else user["experiment"]
    print(user)
    # this is because prior locations for non-negative variables are enterred
    # unlogged:
    if stan_variable.non_negative:
        user["location"] = np.log(user["location"])
    qf = (
        partial(get_lognormal_parameters_from_quantiles, p1=0.01, p2=0.99)
        if stan_variable.non_negative
        else partial(get_normal_parameters_from_quantiles, p1=0.01, p2=0.99)
    )
    user["pct_loc"], user["pct_scale"] = qf(x1=user["pct1"], x2=user["pct99"])
    is_ls = user[["location", "scale"]].notnull().all(axis=1)
    user["used_loc"] = np.where(is_ls, user["location"], user["pct_loc"])
    user["used_scale"] = np.where(is_ls, user["scale"], user["pct_scale"])
    if not user.empty:
        user_loc = user.set_index(["row_id", "col_id"])["used_loc"].unstack()
        user_scale = user.set_index(["row_id", "col_id"])[
            "used_scale"
        ].unstack()
        loc = loc.mask(user_loc.notnull(), user_loc)
        scale = scale.mask(user_scale.notnull(), user_scale)
    return IndPrior2d(stan_variable, loc, scale)


def get_prior_set(upi: UserPriorInput, sv: StanVariableSet) -> PriorSet:
    """Get a PriorSet from a UserPriorInput and StanVariableSet."""
    if upi.dgf_loc is not None and upi.dgf_cov is not None:
        met_order = sv.dgf.ids[0]
        prior_dgf = MultiVariateNormalPrior1d(
            sv.dgf,
            upi.dgf_loc.loc[met_order],
            upi.dgf_cov.loc[met_order, met_order],
        )
    else:
        ip = load_1d_prior(upi.main_table, sv.dgf)
        prior_dgf = MultiVariateNormalPrior1d(
            sv.dgf, ip.location, series_to_diag_df(ip.scale).fillna(0)
        )
    return PriorSet(
        dgf=prior_dgf,
        km=load_1d_prior(upi.main_table, sv.km),
        ki=load_1d_prior(upi.main_table, sv.ki),
        kcat=load_1d_prior(upi.main_table, sv.kcat),
        psi=load_1d_prior(upi.main_table, sv.psi),
        dissociation_constant=load_1d_prior(
            upi.main_table, sv.dissociation_constant
        ),
        transfer_constant=load_1d_prior(upi.main_table, sv.transfer_constant),
        kcat_pme=load_1d_prior(upi.main_table, sv.kcat_pme),
        drain=load_2d_prior(upi.main_table, sv.drain),
        conc_enzyme=load_2d_prior(upi.main_table, sv.conc_enzyme),
        conc_unbalanced=load_2d_prior(upi.main_table, sv.conc_unbalanced),
        conc_pme=load_2d_prior(upi.main_table, sv.conc_pme),
    )
