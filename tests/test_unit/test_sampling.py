"""Unit tests for sampling functions."""

import json
import os

from numpy import ndarray
from numpy.testing import assert_equal

from maud import io, sampling


here = os.path.dirname(__file__)
data_path = os.path.join(here, "..", "data")


def test_get_input_data():
    """Test that the function get_input_data behaves as expected."""
    input_path = os.path.join(data_path, "linear")
    mi = io.load_maud_input_from_toml(input_path)
    with open(os.path.join(input_path, "linear.json"), "r") as f:
        expected_input_data = json.load(f)
    actual_input_data = sampling.get_input_data(mi)
    assert actual_input_data.keys() == expected_input_data.keys()
    for k, v in actual_input_data.items():
        actual = v.tolist() if isinstance(v, ndarray) else v
        expected = expected_input_data[k]
        assert_equal(actual, expected, err_msg=f"{k} different from expected.")
