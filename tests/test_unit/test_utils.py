"""Test functions from the utils module."""


import numpy as np

import maud.utils as utils


def test_get_lognormal_parameters_from_quantiles():
    """Test that the function get_lognormal_parameters_from_quantiles works."""
    inputs = [(0.4, 0.01, 6.8, 0.99)]
    expected = [(0.5003159401539531, 0.608940170915830)]
    for args, (expected_mu, expected_sigma) in zip(inputs, expected):
        mu, sigma = utils.get_lognormal_parameters_from_quantiles(*args)
        assert np.isclose(mu, expected_mu)
        assert np.isclose(sigma, expected_sigma)


def test_get_normal_parameters_from_quantiles():
    """Test that the function get_normal_parameters_from_quantiles works."""
    inputs = [(-3, 0.1, 2, 0.6)]
    expected = [(1.174710655806321, 3.2575440333782133)]
    for args, (expected_mu, expected_sigma) in zip(inputs, expected):
        mu, sigma = utils.get_normal_parameters_from_quantiles(*args)
        assert np.isclose(mu, expected_mu)
        assert np.isclose(sigma, expected_sigma)
