"""Unit tests for sampling functions."""

import json
import os

from numpy.testing import assert_equal

from maud import io, sampling


here = os.path.dirname(__file__)
data_path = os.path.join(here, "..", "data")


def test_get_input_data():
    """Test that the function get_input_data behaves as expected."""
    input_path = os.path.join(data_path, "linear")
    mi = io.load_maud_input_from_toml(input_path)
    expected = json.load(open(os.path.join(input_path, "linear.json"), "r"))
    actual = sampling.get_input_data(mi)
    assert actual.keys() == expected.keys()
    for k in actual.keys():
        assert_equal(actual[k], expected[k], err_msg=f"{k} is different from expected.")
