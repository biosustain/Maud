"""Unit tests for io functions."""

import os

import maud.io as io


here = os.path.dirname(__file__)
data_path = os.path.join(here, "../data")


def test_load_maud_input_from_toml():
    """Test that the function load_maud_input_from_toml behaves as expected."""
    expected_stan_codes = {
        "metabolite": {"M1": 1, "M2": 2},
        "metabolite_in_compartment": {"M1_e": 1, "M2_e": 2, "M1_c": 3, "M2_c": 4},
        "balanced_mic": {"M1_c": 3, "M2_c": 4},
        "unbalanced_mic": {"M1_e": 1, "M2_e": 2},
        "reaction": {"r1": 2, "r2": 1, "r3": 3},
        "experiment": {"condition_1": 1, "condition_2": 2},
        "enzyme": {"r1": 2, "r2": 1, "r3": 3},
        "drain": None,
    }
    mi = io.load_maud_input_from_toml(os.path.join(data_path, "linear.toml"))
    assert mi.kinetic_model.reactions["r1"].stoichiometry == {"M1_e": -1, "M1_c": 1}
    assert "r1" in map(lambda p: p.enzyme_id, mi.priors["kcats"])
    assert (
        mi.experiments["condition_1"].measurements["metabolite"]["M1_c"].target_type
        == "metabolite"
    )
    assert mi.stan_codes == expected_stan_codes
