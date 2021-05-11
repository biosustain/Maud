"""Unit tests for io functions."""

import os

import maud.io as io


here = os.path.dirname(__file__)
data_path = os.path.join(here, "../data")


def test_load_maud_input_from_toml():
    """Test that the function load_maud_input_from_toml behaves as expected."""
    expected_stan_coords = {
        "metabolites": ["M1", "M2"],
        "mics": ["M1_c", "M1_e", "M2_c", "M2_e"],
        "balanced_mics": ["M1_c", "M2_c"],
        "unbalanced_mics": ["M1_e", "M2_e"],
        "reactions": ["r1", "r2", "r3"],
        "experiments": ["condition_1", "condition_2"],
        "enzymes": ["r1", "r2", "r3"],
        "phos_enzs": [],
        "drains": [],
    }
    mi = io.load_maud_input_from_toml(os.path.join(data_path, "linear"))
    r1 = [r for r in mi.kinetic_model.reactions if r.id == "r1"][0]
    assert r1.stoichiometry == {"M1_e": -1, "M1_c": 1}
    assert "r1" in mi.priors.priors_kcat.location.reset_index()["enzyme_id"].values
    for key, expected_codes in expected_stan_coords.items():
        assert mi.stan_coords.__dict__[key] == expected_codes
