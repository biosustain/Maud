"""Unit tests for io functions."""

import json

import importlib_resources
import numpy as np
from numpy.testing import assert_equal

from maud.data.example_inputs import linear, linear_multidgf
from maud.loading_maud_inputs import load_maud_input


def test_load_maud_input():
    """Test that the function load_maud_input behaves as expected."""
    expected_var_ids = {
        "dgf": [["M1", "M2"]],
        "conc_unbalanced_train": [
            ["condition1", "condition2"],
            ["M1_e", "M2_e"],
        ],
        "conc_pme_train": [["condition1", "condition2"], []],
        "dissociation_constant": [["r1_M2_c_activation", "r2_M1_c_inhibition"]],
    }
    linear_files = importlib_resources.files(linear)
    mi = load_maud_input(data_path=str(linear_files))  # path 0 is package
    r1 = next(r for r in mi.kinetic_model.reactions if r.id == "r1")
    assert r1.stoichiometry == {"M1_e": -1, "M1_c": 1}
    assert "r1_r1" in mi.parameters.kcat.ids[0]
    for var_name, expected in expected_var_ids.items():
        assert getattr(mi.parameters, var_name).ids == expected
    actual_stan_input = mi.stan_input_train
    with linear_files.joinpath("expected_stan_input.json").open("r") as f:
        expected_stan_input = json.load(f)
    assert set(actual_stan_input.keys()) == set(expected_stan_input.keys())
    for k, v in actual_stan_input.items():
        actual = v.tolist() if isinstance(v, np.ndarray) else v
        expected = expected_stan_input[k]
        assert_equal(actual, expected, err_msg=f"{k} different from expected.")


def test_load_multidgf():
    """Test that the multidgf input loads correctly."""
    files = importlib_resources.files(linear_multidgf)
    mi = load_maud_input(data_path=str(files))  # path 0 is package
    assert mi.inits_dict["dgf_free"] == [-10.0, -32.0]
