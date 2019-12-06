"""Unit tests for sampling functions."""

import os

from numpy.testing import assert_equal

from maud import io, sampling


here = os.path.dirname(__file__)
data_path = os.path.join(here, "../data")


def test_get_input_data():
    """Test that the function get_input_data behaves as expected."""
    expected = {
        "N_balanced": 2,
        "N_unbalanced": 2,
        "N_kinetic_parameters": 11,
        "N_reaction": 3,
        "N_enzyme": 3,
        "N_experiment": 2,
        "N_flux_measurement": 2,
        "N_conc_measurement": 4,
        "stoichiometric_matrix": [[-1, 0, 0], [0, 0, 1], [1, -1, 0], [0, 1, -1]],
        "experiment_yconc": [1, 1, 2, 2],
        "metabolite_yconc": [3, 4, 3, 4],
        "yconc": [0.8, 1.5, 0.7, 1.4],
        "sigma_conc": [0.1, 0.1, 0.1, 0.1],
        "experiment_yflux": [1, 2],
        "reaction_yflux": [3, 3],
        "yflux": [0.29, 0.21],
        "sigma_flux": [0.1, 0.1],
        "prior_loc_formation_energy": [-1.0, -2.0, -1.0, -2.0],
        "prior_scale_formation_energy": [0.05, 0.05, 0.05, 0.05],
        "prior_loc_kinetic_parameters": [
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
        ],
        "prior_scale_kinetic_parameters": [
            0.6,
            0.6,
            0.6,
            0.6,
            0.6,
            0.6,
            0.6,
            0.6,
            0.6,
            0.6,
            0.6,
        ],
        "prior_loc_unbalanced": [[2.0, 2.0], [1.0, 1.0]],
        "prior_scale_unbalanced": [[0.05, 0.05], [0.05, 0.05]],
        "prior_loc_enzyme": [[1.0, 1.5], [1.0, 1.5], [1.0, 1.5]],
        "prior_scale_enzyme": [[0.05, 0.05], [0.05, 0.05], [0.05, 0.05]],
        "as_guess": [0.1, 0.1],
        "rtol": 1e-09,
        "ftol": 1e-06,
        "steps": 1000000000,
        "LIKELIHOOD": 1,
    }
    toml_input_path = os.path.join(data_path, "linear.toml")
    mi = io.load_maud_input_from_toml(toml_input_path)
    actual, _ = sampling.get_input_data(mi, 1e-06, 1e-09, int(1e9), 1)
    assert actual.keys() == expected.keys()
    for k in actual.keys():
        assert_equal(actual[k], expected[k])
