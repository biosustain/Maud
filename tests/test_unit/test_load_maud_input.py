"""Unit tests for io functions."""

import json
import os

import numpy as np
from numpy.testing import assert_equal

from maud.loading_maud_inputs import load_maud_input


HERE = os.path.dirname(__file__)
LINEAR_PATH = os.path.join(HERE, "..", "data", "linear")
EXPECTED_STAN_INPUT_PATH = os.path.join(
    LINEAR_PATH, "expected_stan_input.json"
)


def test_load_maud_input():
    """Test that the function load_maud_input behaves as expected."""
    expected_var_ids = {
        "dgf": [["M1", "M2"]],
        "conc_unbalanced": [["condition1", "condition2"], ["M1_e", "M2_e"]],
        "conc_phos": [["condition1", "condition2"], []],
        "dissociation_constant": [["r1_M2_c", "r2_M1_c"]],
    }
    mi = load_maud_input(data_path=LINEAR_PATH)
    r1 = next(r for r in mi.kinetic_model.reactions if r.id == "r1")
    assert r1.stoichiometry == {"M1_e": -1, "M1_c": 1}
    assert "r1_r1" in mi.priors.kcat.location.index
    for var_name, expected in expected_var_ids.items():
        assert var_name in mi.stan_variable_set.__dataclass_fields__
        assert getattr(mi.stan_variable_set, var_name).ids == expected
    actual_stan_input = mi.stan_input_train.stan_input_dict
    with open(EXPECTED_STAN_INPUT_PATH, "r") as f:
        expected_stan_input = json.load(f)
    assert set(actual_stan_input.keys()) == set(expected_stan_input.keys())
    for k, v in actual_stan_input.items():
        actual = v.tolist() if isinstance(v, np.ndarray) else v
        expected = expected_stan_input[k]
        assert_equal(actual, expected, err_msg=f"{k} different from expected.")
