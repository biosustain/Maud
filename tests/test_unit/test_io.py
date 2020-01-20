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
        "kinetic_parameter": {
            "r1_Ka": 1,
            "r1_Kp": 2,
            "r1_Kcat1": 3,
            "r1_dissociation_constant_r_M2_c": 4,
            "r1_transfer_constant": 5,
            "r2_Ka": 6,
            "r2_Kp": 7,
            "r2_Kcat1": 8,
            "r2_dissociation_constant_t_M1_c": 9,
            "r2_transfer_constant": 10,
            "r2_inhibition_constant_M2_c": 11,
            "r3_Ka": 12,
            "r3_Kp": 13,
            "r3_Kcat1": 14,
        },
        "reaction": {"r1": 1, "r2": 2, "r3": 3},
        "experiment": {"condition_1": 1, "condition_2": 2},
        "enzyme": {"r1": 1, "r2": 2, "r3": 3},
    }
    mi = io.load_maud_input_from_toml(os.path.join(data_path, "linear.toml"))
    assert mi.kinetic_model.reactions["r1"].stoichiometry == {"M1_e": -1, "M1_c": 1}
    print(mi.priors["r1_Kcat1"].target_id)
    assert mi.priors["r1_Kcat1"].target_id == "Kcat1"
    assert (
        mi.experiments["condition_1"].measurements["metabolite"]["M1_c"].target_type
        == "metabolite"
    )
    assert mi.stan_codes == expected_stan_codes
