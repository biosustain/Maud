"""General purpose utility functions."""

from typing import Any, Dict, Hashable, List

import numpy as np
import pandas as pd
import sympy as sp
from depinfo import print_dependencies
from scipy.stats import norm


def join_str_cols(df: pd.DataFrame, sep: str, name=None) -> pd.Series:
    """Join the columns of a dataframe into a Series separated by sep."""
    out = df.apply(sep.join, axis=1).rename(name)
    assert isinstance(out, pd.Series)
    return out


def load_df(read_csv_input, **kwargs) -> pd.DataFrame:
    """Wrap pd.read_csv, ensuring that a dataframe is returned."""
    return check_is_df(pd.read_csv(read_csv_input, **kwargs))


def show_versions():
    """Print dependency information."""
    print_dependencies("maud")


def recursively_flatten_list(o: List) -> List:
    """Recursively flatten a nested list."""
    gather = []
    for item in o:
        if isinstance(item, List):
            gather.extend(recursively_flatten_list(item))
        else:
            gather.append(item)
    return gather


def check_is_df(maybe_df: Any) -> pd.DataFrame:
    """Assert that an object is a pandas DataFrame, then return it.

    This is useful for keeping mypy happy.
    """
    assert isinstance(maybe_df, pd.DataFrame)
    return maybe_df


def series_to_diag_df(s: pd.Series) -> pd.DataFrame:
    """Turn a Series into a Dataframe with it on the diagonal.

    The off-diagonal cells are nan.
    """
    out = pd.DataFrame(np.nan, index=s.index, columns=s.index)
    np.fill_diagonal(out.values, s.values)
    return out


def read_with_fallback(k: Hashable, d: Dict, default: Any):
    """Get item k from d if it is available, otherwise return default."""
    return d[k] if k in d.keys() else default


def codify(lx: List) -> Dict[str, int]:
    """Turn a list of strings into a dictionary mapping them to integers."""
    return dict(zip(lx, range(1, len(lx) + 1)))


def get_lognormal_parameters_from_quantiles(x1, p1, x2, p2):
    """Find parameters for a lognormal distribution from two quantiles.

    i.e. get mu and sigma such that if X ~ lognormal(mu, sigma), then pr(X <
    x1) = p1 and pr(X < x2) = p2.

    """
    logx1 = np.log(x1)
    logx2 = np.log(x2)
    denom = norm.ppf(p2) - norm.ppf(p1)
    sigma = (logx2 - logx1) / denom
    mu = (logx1 * norm.ppf(p2) - logx2 * norm.ppf(p1)) / denom
    return mu, sigma


def get_normal_parameters_from_quantiles(x1, p1, x2, p2):
    """Find parameters for a normal distribution from two quantiles.

    i.e. get mu and sigma such that if X ~ normal(mu, sigma), then pr(X <
    x1) = p1 and pr(X < x2) = p2.

    """
    denom = norm.ppf(p2) - norm.ppf(p1)
    sigma = (x2 - x1) / denom
    mu = (x1 * norm.ppf(p2) - x2 * norm.ppf(p1)) / denom
    return mu, sigma


def get_null_space(a, rtol=1e-5):
    """Calulate the null space of a matrix."""
    u, s, v = np.linalg.svd(a)
    rank = (s > rtol * s[0]).sum()
    return v[rank:].T.copy()


def get_left_nullspace(matrix: pd.DataFrame, atol=1e-13, rtol=0.0):
    """Compute an approximate basis for the null space (kernel) of a matrix.

    The algorithm used by this function is based on the singular value
    decomposition of the given matrix.

    Parameters
    ----------
    matrix : ndarray
        The matrix should be at most 2-D.  A 1-D array with length k
        will be treated as a 2-D with shape (1, k)
    atol : float
        The absolute tolerance for a zero singular value.  Singular values
        smaller than ``atol`` are considered to be zero.
    rtol : float
        The relative tolerance for a zero singular value.  Singular values less
        than the relative tolerance times the largest singular value are
        considered to be zero.

    Notes
    -----
    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
        tol = max(atol, rtol * smax)
    Singular values smaller than ``tol`` are considered to be zero.

    Returns
    -------
    ndarray
        If ``matrix`` is an array with shape (m, k), then the returned
        nullspace will be an array with shape ``(k, n)``, where n is the
        estimated dimension of the nullspace.

    References
    ----------
    Adapted from:
    https://scipy.github.io/old-wiki/pages/Cookbook/RankNullspace.html
    and then taken from from
    https://github.com/opencobra/memote/blob/develop/src/memote/support/consistency_helpers.py#L163

    """  # noqa: D402
    matrix = np.atleast_2d(matrix)
    _, sigma, vh = np.linalg.svd(matrix.T)
    tol = max(atol, rtol * sigma[0])
    num_nonzero = (sigma >= tol).sum()
    return vh[num_nonzero:].conj().T


def get_rref(mat):
    """Return reduced row echelon form of a matrix."""
    return sp.Matrix(mat).rref(iszerofunc=lambda x: abs(x) < 1e-10)[0]
