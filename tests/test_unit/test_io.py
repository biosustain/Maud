"""Unit tests for io functions."""

import os

import maud.io as io


here = os.path.dirname(__file__)
data_path = os.path.join(here, "../data")


def test_load_maud_input_from_toml():
    """Test that the function load_maud_input_from_toml behaves as expected."""
    expected_stan_codes = {
        "metabolite_codes": {"M1": 1, "M2": 2},
        "mic_codes": {"M1_e": 1, "M2_e": 2, "M1_c": 3, "M2_c": 4},
        "balanced_mic_codes": {"M1_c": 3, "M2_c": 4},
        "unbalanced_mic_codes": {"M1_e": 1, "M2_e": 2},
        "reaction_codes": {"r1": 1, "r2": 2, "r3": 3},
        "experiment_codes": {"condition_1": 1, "condition_2": 2},
        "enzyme_codes": {"r1": 1, "r2": 2, "r3": 3},
        "drain_codes": {},
    }
    mi = io.load_maud_input_from_toml(os.path.join(data_path, "linear.toml"))
    assert mi.kinetic_model.reactions["r1"].stoichiometry == {"M1_e": -1, "M1_c": 1}
    assert "r1" in map(lambda p: p.enzyme_id, mi.priors.kcat_priors)
    exp = [e for e in mi.experiments.experiments if e.id == "condition_1"][0]
    assert exp.measurements["metabolite"]["M1_c"].target_type == "metabolite"
    assert mi.stan_codes.__dict__ == expected_stan_codes
