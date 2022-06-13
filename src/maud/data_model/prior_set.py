"""Definitions of prior-related things.

The most important is PriorSet, which defines the set of priors that any
MaudInput must have.

"""

from typing import Optional

import numpy as np
import pandas as pd
from numpy.linalg import LinAlgError
from pydantic.class_validators import root_validator, validator
from pydantic.dataclasses import dataclass

from maud.data_model.stan_variable_set import StanVariable


class PriorConfig:
    """Config allowing priors to contain pandas objects."""

    arbitrary_types_allowed = True


@dataclass(config=PriorConfig)
class UserPriorInput:
    """The user's prior input, consisting of a dataframe and maybe dgf info."""

    main_table: pd.DataFrame
    dgf_loc: Optional[pd.Series]
    dgf_cov: Optional[pd.DataFrame]


@dataclass(config=PriorConfig)
class IndPrior1d:
    """Independent location/scale prior for a 1D parameter."""

    stan_variable: StanVariable
    location: pd.Series
    scale: pd.Series

    @root_validator
    def loc_and_scale_indexes_must_match(cls, values):
        """Check that location and scale have the same index."""
        assert values["location"].index.equals(
            values["scale"].index
        ), "Location index doesn't match scale index."
        return values


@dataclass(config=PriorConfig)
class IndPrior2d:
    """Independent location/scale prior for a 2D parameter."""

    stan_variable: StanVariable
    location: pd.DataFrame
    scale: pd.DataFrame

    @root_validator
    def loc_and_scale_indexes_and_columns_must_match(cls, values):
        """Check that location and scale have the same index and columns."""
        assert values["location"].index.equals(
            values["scale"].index
        ), "Location index doesn't match scale index."
        assert values["location"].columns.equals(
            values["scale"].columns
        ), "Location columns doen't match scale columns."
        return values


@dataclass(config=PriorConfig)
class MultiVariateNormalPrior1d:
    """Prior Location vector and covariance matrix for a 1D parameter."""

    stan_variable: StanVariable
    location: pd.Series
    covariance_matrix: pd.DataFrame

    @root_validator
    def location_index_must_agree_with_stan_variable(cls, values):
        """Check that location index is the same as stan variable ids."""
        location_index = values["location"].index.tolist()
        sv_index = values["stan_variable"].ids[0]
        for ix_loc, ix_sv in zip(location_index, sv_index):
            assert ix_loc == ix_sv
        return values

    @root_validator
    def cov_and_loc_indexes_must_match(cls, values):
        """Check that location and cov agree."""
        assert values["location"].index.equals(
            values["covariance_matrix"].index
        ), "Location index doesn't match covariance matrix index."
        assert values["location"].index.equals(
            values["covariance_matrix"].columns
        ), "Location columns don't match covariance matrix columns."
        return values

    @validator("covariance_matrix")
    def cov_matrix_must_be_pos_def(cls, v):
        """Check that covariance matrix is positive definite."""
        try:
            np.linalg.cholesky(v.values)
        except LinAlgError as e:
            raise ValueError(
                "Covariance matrix is not positive definite"
            ) from e
        return v


@dataclass
class PriorSet:
    """Object containing all priors for a MaudInput."""

    dgf: MultiVariateNormalPrior1d
    km: IndPrior1d
    kcat: IndPrior1d
    kcat_pme: IndPrior1d
    ki: IndPrior1d
    psi: IndPrior1d
    dissociation_constant: IndPrior1d
    transfer_constant: IndPrior1d
    conc_unbalanced: IndPrior2d
    drain: IndPrior2d
    conc_enzyme: IndPrior2d
    conc_pme: IndPrior2d

    @root_validator
    def missing_priors_are_not_allowed(cls, values):
        """Check that there are no missing priors."""
        for k, v in values.items():
            if isinstance(v.location, pd.DataFrame):
                missing_loc = v.location.loc[:, v.location.isnull().any(0)]
            else:
                missing_loc = v.location[v.location.isnull()]
            assert (
                missing_loc.empty
            ), f"Missing {v.location.isnull().sum().sum()} {k}: {missing_loc}"
        return values
