"""Definitions of priors."""

import math
from typing import List, Optional, Union

import numpy as np
import pandas as pd
from numpy.linalg import LinAlgError
from pydantic import (
    BaseModel,
    NonNegativeFloat,
    PositiveFloat,
    computed_field,
    field_validator,
    model_validator,
)

from maud.data_model.hardcoding import ID_SEPARATOR
from maud.data_model.id_component import IdComponent
from maud.data_model.parameter_input import (
    ParameterInputAtom,
    ParameterInputMVN,
)
from maud.utility_functions import (
    get_lognormal_parameters_from_quantiles,
    get_normal_parameters_from_quantiles,
)


def get_pia_loc(pia: ParameterInputAtom, non_negative: bool) -> float:
    """Get the location from a parameter input atom."""
    if pia.location is not None:
        return pia.location
    elif pia.exploc is not None:
        return np.log(pia.exploc)
    elif pia.pct1 is not None:
        f = (
            get_lognormal_parameters_from_quantiles
            if non_negative
            else get_normal_parameters_from_quantiles
        )
        loc, _ = f(pia.pct1, 0.01, pia.pct99, 0.99)
        return loc
    else:
        raise ValueError(f"Incorrectly specified Prior input atom: {pia}")


def get_pia_scale(pia: ParameterInputAtom, non_negative: bool) -> float:
    """Get the scale from a parameter input atom."""
    if pia.scale is not None:
        return pia.scale
    elif pia.pct1 is not None:
        f = (
            get_lognormal_parameters_from_quantiles
            if non_negative
            else get_normal_parameters_from_quantiles
        )
        _, scale = f(pia.pct1, 0.01, pia.pct99, 0.99)
        return scale
    else:
        raise ValueError(f"Incorrectly specified Prior input atom: {pia}")


class IndPrior1d(BaseModel):
    """Independent location/scale prior for a 1D parameter."""

    user_input: Optional[List[ParameterInputAtom]]
    ids: List[List[str]]
    id_components: List[List[IdComponent]]
    non_negative: bool
    default_loc: float
    default_scale: PositiveFloat

    @computed_field
    def location(self) -> List[float]:
        """Add the location field."""
        if len(self.ids[0]) == 0:
            return []
        loc_series = pd.Series(self.default_loc, index=self.ids)
        if self.user_input is not None:
            for pia in self.user_input:
                loc_i = get_pia_loc(pia, non_negative=self.non_negative)
                ids_i = [
                    ID_SEPARATOR.join([getattr(pia, c) for c in idci])
                    for idci in self.id_components
                ]
                if ids_i[0] in loc_series.index:
                    loc_series.loc[ids_i[0]] = loc_i
        return loc_series.tolist()

    @computed_field
    def scale(self) -> List[PositiveFloat]:
        """Add the scale field."""
        if len(self.ids[0]) == 0:
            return []
        scale_series = pd.Series(self.default_scale, index=self.ids)
        if self.user_input is not None:
            for pia in self.user_input:
                scale_i = get_pia_scale(pia, non_negative=self.non_negative)
                ids_i = [
                    ID_SEPARATOR.join([getattr(pia, c) for c in idci])
                    for idci in self.id_components
                ]
                if ids_i[0] in scale_series.index:
                    scale_series.loc[ids_i[0]] = scale_i
        return scale_series.tolist()

    @field_validator("location")
    def no_null_locations(cls, v):
        """Check that locations are non-null."""
        for x in v:
            if math.isnan(x):
                raise ValueError("Location cannot be null.")
        return v

    @field_validator("scale")
    def no_null_scales(cls, v):
        """Check that scales are non-null."""
        for x in v:
            if math.isnan(x):
                raise ValueError("Scale cannot be null.")
        return v

    @field_validator("scale")
    def scales_are_positive(cls, v):
        """Check that scales are positive."""
        if len(v) > 0 and any(x <= 0 for x in v):
            raise ValueError("Scale parameter must be positive.")
        return v

    @model_validator(mode="after")
    def lengths_match(self):
        """Check that location and scale have the same index."""
        n_locs = len(self.location)
        n_scales = len(self.scale)
        if n_locs != n_scales:
            raise ValueError(
                "Location, scale and ids must have the same length."
            )
        return self


class IndPrior2d(BaseModel):
    """Independent location/scale prior for a 2D parameter."""

    user_input: Optional[List[ParameterInputAtom]]
    ids: List[List[str]]
    id_components: List[List[IdComponent]]
    non_negative: bool
    default_loc: float
    default_scale: PositiveFloat

    @computed_field
    def location(self) -> List[List[float]]:
        """Add the location field."""
        if any(len(ids_i) == 0 for ids_i in self.ids):
            return [[]]
        loc_df = pd.DataFrame(
            self.default_loc, index=self.ids[0], columns=self.ids[1]
        )
        if self.user_input is not None:
            for pia in self.user_input:
                loc_i = get_pia_loc(pia, non_negative=self.non_negative)
                ids_i = [
                    ID_SEPARATOR.join([getattr(pia, c) for c in idci])
                    for idci in self.id_components
                ]
                if ids_i[0] in loc_df.index and ids_i[1] in loc_df.columns:
                    loc_df.loc[ids_i[0], ids_i[1]] = loc_i
        return loc_df.values.tolist()

    @computed_field
    def scale(self) -> List[List[PositiveFloat]]:
        """Add the scale field."""
        if any(len(ids_i) == 0 for ids_i in self.ids):
            return [[]]
        scale_df = pd.DataFrame(
            self.default_scale, index=self.ids[0], columns=self.ids[1]
        )
        if self.user_input is not None:
            for pia in self.user_input:
                scale_i = get_pia_scale(pia, non_negative=self.non_negative)
                ids_i = [
                    ID_SEPARATOR.join([getattr(pia, c) for c in idci])
                    for idci in self.id_components
                ]
                if ids_i[0] in scale_df.index and ids_i[1] in scale_df.columns:
                    scale_df.loc[ids_i[0], ids_i[1]] = scale_i
        return scale_df.values.tolist()

    @field_validator("location")
    def no_null_locations(cls, v):
        """Check that locations are non-null."""
        for vi in v:
            for x in vi:
                if math.isnan(x):
                    raise ValueError("Location cannot be null.")
        return v

    @field_validator("scale")
    def no_null_scales(cls, v):
        """Check that scales are non-null."""
        for vi in v:
            for x in vi:
                if math.isnan(x):
                    raise ValueError("Scale cannot be null.")
        return v

    @field_validator("scale")
    def scales_are_positive(cls, v):
        """Check that scales are all positive."""
        for s in v:
            if len(s) > 0 and any(x <= 0 for x in s):
                raise ValueError("Scale parameter must be positive.")
        return v

    @model_validator(mode="after")
    def lengths_are_correct(self):
        """Check that ids, location and scale have correct length."""
        loc = self.location
        scale = self.scale
        if not len(loc) == len(scale):
            raise ValueError("First dimension length incorrect.")
        for i, (loc_i, scale_i) in enumerate(zip(loc, scale)):
            if not len(loc_i) == len(scale_i):
                raise ValueError(f"Length mismatch at index {i}.")
        return self


class PriorMVN(BaseModel):
    """Prior Location vector and covariance matrix for a 1D parameter."""

    user_input: Optional[Union[List[ParameterInputAtom], ParameterInputMVN]]
    ids: List[List[str]]
    id_components: List[List[IdComponent]]
    non_negative: bool
    default_loc: float
    default_scale: PositiveFloat

    @computed_field
    def location(self) -> List[float]:
        """Add the location field."""
        ids = self.ids[0]
        loc_series = pd.Series(self.default_loc, index=self.ids)
        if isinstance(self.user_input, ParameterInputMVN):
            loc_series = pd.Series(
                self.user_input.mean_vector, index=self.user_input.ids
            ).reindex(ids)
        elif self.user_input is not None:
            for pia in self.user_input:
                ids_i = [
                    ID_SEPARATOR.join([getattr(pia, c) for c in idci])
                    for idci in self.id_components
                ]
                loc_i = get_pia_loc(pia, non_negative=self.non_negative)
                loc_series.loc[ids_i[0]] = loc_i
        return loc_series.tolist()

    @computed_field
    def covariance_matrix(self) -> List[List[NonNegativeFloat]]:
        """Add the covariance_matrix field."""
        ids = self.ids[0]
        cov_df = pd.DataFrame(
            np.diagflat(np.tile(self.default_scale, len(ids))),
            index=ids,
            columns=ids,
        )
        if isinstance(self.user_input, ParameterInputMVN):
            cov_df = (
                pd.DataFrame(
                    self.user_input.covariance_matrix,
                    index=self.user_input.ids,
                    columns=self.user_input.ids,
                )
                .reindex(ids)
                .reindex(columns=ids)
            )
        elif self.user_input is not None:
            for pia in self.user_input:
                ids_i = [
                    ID_SEPARATOR.join([getattr(pia, c) for c in idci])
                    for idci in self.id_components
                ]
                cov_ii = get_pia_scale(pia, non_negative=self.non_negative)
                cov_df.loc[ids_i[0], ids_i[0]] = cov_ii
        return cov_df.values.tolist()

    @field_validator("covariance_matrix")
    def no_null_covariances(cls, v):
        """Check that scales are non-null."""
        for vi in v:
            for x in vi:
                if math.isnan(x):
                    raise ValueError("Covariance cannot be null.")
        return v

    @field_validator("location")
    def no_null_locations(cls, v):
        """Check that no locations are nans."""
        for x in v:
            if math.isnan(x):
                raise ValueError("Location cannot be null.")
        return v

    @field_validator("covariance_matrix")
    def cov_matrix_is_pos_def(cls, v):
        """Check that covariance matrix is positive definite."""
        try:
            np.linalg.cholesky(v)
        except LinAlgError as e:
            raise ValueError(
                "Covariance matrix is not positive definite"
            ) from e
        return v

    @model_validator(mode="after")
    def lengths_are_correct(self):
        """Check that ids, location and cov matrix have correct lengths."""
        loc = self.location
        cov = self.covariance_matrix
        if not all(len(x) == len(cov[0]) for x in cov):
            raise ValueError(
                "All elements of covariance matrix must have same length."
            )
        if not len(loc) == len(cov[0]):
            raise ValueError("First dimension length incorrect.")
        return self


Prior = Union[IndPrior1d, IndPrior2d, PriorMVN]
