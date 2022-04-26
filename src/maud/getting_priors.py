"""Functions for creating *Prior* and PriorSet objects."""

import pandas as pd

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
    loc = pd.Series(stan_variable.default_loc, index=stan_variable.ids[0])
    scale = pd.Series(stan_variable.default_scale, index=stan_variable.ids[0])
    user_df_for_param = user_df.loc[
        lambda df: df["parameter"] == stan_variable.name
    ]
    for _, row in user_df_for_param.iterrows():
        if row[["location", "scale"]].notnull().all():
            loc.loc[row["row_id"]] = row["location"]
            scale.loc[row["row_id"]] = row["scale"]
        elif row["pct1"].notnull() and row["pct99"].notnull():
            qf = (
                get_lognormal_parameters_from_quantiles
                if stan_variable.non_negative
                else get_normal_parameters_from_quantiles
            )
            loc.loc[row["row_id"]], scale.loc[row["row_id"]] = qf(
                row["pct1"], 0.01, row["pct99"], 0.99
            )
    return IndPrior1d(stan_variable, loc, scale)


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
    user_df_for_param = user_df.loc[
        lambda df: df["parameter"] == stan_variable.name
    ]
    for _, row in user_df_for_param.iterrows():
        if row[["location", "scale"]].notnull().all():
            loc.loc[row["row_id"], row["col_id"]] = row["location"]
            scale.loc[row["row_id"], row["col_id"]] = row["scale"]
        elif row["pct1"].notnull() and row["pct99"].notnull():
            qf = (
                get_lognormal_parameters_from_quantiles
                if stan_variable.non_negative
                else get_normal_parameters_from_quantiles
            )
            (
                loc.loc[row["row_id"], row["col_id"]],
                scale.loc[row["row_id"], row["col_id"]],
            ) = qf(row["pct1"], 0.01, row["pct99"], 0.99)
    return IndPrior2d(stan_variable, loc, scale)


def get_prior_set(upi: UserPriorInput, sv: StanVariableSet) -> PriorSet:
    """Get a PriorSet from a UserPriorInput and StanVariableSet."""
    if upi.dgf_loc is not None and upi.dgf_cov is not None:
        prior_dgf = MultiVariateNormalPrior1d(sv.dgf, upi.dgf_loc, upi.dgf_cov)
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
        kcat_phos=load_1d_prior(upi.main_table, sv.kcat_phos),
        drain=load_2d_prior(upi.main_table, sv.drain),
        conc_enzyme=load_2d_prior(upi.main_table, sv.conc_enzyme),
        conc_unbalanced=load_2d_prior(upi.main_table, sv.conc_unbalanced),
        conc_phos=load_2d_prior(upi.main_table, sv.conc_phos),
    )
