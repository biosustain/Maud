"""Definitions of priors."""

import math
from typing import List, Union

import numpy as np
from numpy.linalg import LinAlgError
from pydantic.class_validators import root_validator, validator
from pydantic.dataclasses import dataclass


@dataclass
class IndPriorAtom:
    """Prior for a single quantity."""

    location: float
    scale: float


@dataclass
class IndPrior1d:
    """Independent location/scale prior for a 1D parameter."""

    location: List[float]
    scale: List[float]

    @validator("location")
    def no_null_locations(cls, v):
        """Check that locations are non-null."""
        for x in v:
            if math.isnan(x):
                raise ValueError("Location cannot be null.")
        return v

    @validator("scale")
    def no_null_scales(cls, v):
        """Check that scales are non-null."""
        for x in v:
            if math.isnan(x):
                raise ValueError("Scale cannot be null.")
        return v

    @validator("scale")
    def scales_are_positive(cls, v):
        """Check that scales are positive."""
        if len(v) > 0 and any(x <= 0 for x in v):
            raise ValueError("Scale parameter must be positive.")
        return v

    @root_validator(pre=False, skip_on_failure=True)
    def lengths_match(cls, values):
        """Check that location and scale have the same index."""
        n_locs = len(values["location"])
        n_scales = len(values["scale"])
        if n_locs != n_scales:
            raise ValueError(
                "Location, scale and ids must have the same length."
            )
        return values


@dataclass
class IndPrior2d:
    """Independent location/scale prior for a 2D parameter."""

    location: List[List[float]]
    scale: List[List[float]]

    @validator("location")
    def no_null_locations(cls, v):
        """Check that locations are non-null."""
        for vi in v:
            for x in vi:
                if math.isnan(x):
                    raise ValueError("Location cannot be null.")
        return v

    @validator("scale")
    def no_null_scales(cls, v):
        """Check that scales are non-null."""
        for vi in v:
            for x in vi:
                if math.isnan(x):
                    raise ValueError("Scale cannot be null.")
        return v

    @validator("scale")
    def scales_are_positive(cls, v):
        """Check that scales are all positive."""
        for s in v:
            if len(s) > 0 and any(x <= 0 for x in s):
                raise ValueError("Scale parameter must be positive.")
        return v

    @root_validator(pre=False, skip_on_failure=True)
    def lengths_are_correct(cls, values):
        """Check that ids, location and scale have correct length."""
        loc = values["location"]
        scale = values["scale"]
        if not len(loc) == len(scale):
            raise ValueError("First dimension length incorrect.")
        for i, (loc_i, scale_i) in enumerate(zip(loc, scale)):
            if not len(loc_i) == len(scale_i):
                raise ValueError(f"Length mismatch at index {i}.")
        return values


@dataclass
class PriorMVN:
    """Prior Location vector and covariance matrix for a 1D parameter."""

    location: List[float]
    covariance_matrix: List[List[float]]

    @validator("covariance_matrix")
    def no_null_covariances(cls, v):
        """Check that scales are non-null."""
        for vi in v:
            for x in vi:
                if math.isnan(x):
                    raise ValueError("Covariance cannot be null.")
        return v

    @validator("location")
    def no_null_locations(cls, v):
        """Check that no locations are nans."""
        for x in v:
            if math.isnan(x):
                raise ValueError("Location cannot be null.")
        return v

    @validator("covariance_matrix")
    def cov_matrix_is_pos_def(cls, v):
        """Check that covariance matrix is positive definite."""
        try:
            np.linalg.cholesky(v)
        except LinAlgError as e:
            raise ValueError(
                "Covariance matrix is not positive definite"
            ) from e
        return v

    @root_validator(pre=False, skip_on_failure=True)
    def lengths_are_correct(cls, values):
        """Check that ids, location and cov matrix have correct lengths."""
        loc = values["location"]
        cov = values["covariance_matrix"]
        if not all(len(x) == len(cov[0]) for x in cov):
            raise ValueError(
                "All elements of covariance matrix must have same length."
            )
        if not len(loc) == len(cov[0]):
            raise ValueError("First dimension length incorrect.")
        return values


Prior = Union[IndPrior1d, IndPrior2d, PriorMVN]
