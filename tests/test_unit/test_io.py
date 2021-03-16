"""Unit tests for io functions."""

import os

import maud.io as io


here = os.path.dirname(__file__)
data_path = os.path.join(here, "../data")


def test_load_maud_input_from_toml():
    """Test that the function load_maud_input_from_toml behaves as expected."""
    expected_stan_codes = {
        "metabolite_codes": ["M1", "M2"],
        "mic_codes": ["M1_e", "M2_e", "M1_c", "M2_c"],
        "balanced_mic_codes": ["M1_c", "M2_c"],
        "unbalanced_mic_codes": ["M1_e", "M2_e"],
        "reaction_codes": ["r1", "r2", "r3"],
        "experiment_codes": ["condition_1", "condition_2"],
        "enzyme_codes": ["r1", "r2", "r3"],
        "phos_enz_codes": [],
        "drain_codes": [],
    }
    mi = io.load_maud_input_from_toml(os.path.join(data_path, "linear"))
    r1 = [r for r in mi.kinetic_model.reactions if r.id == "r1"][0]
    experiment_ids = [m.experiment_id for m in mi.measurements]
    expected_experiment_ids = ["condition_1", "condition_2"]
    assert r1.stoichiometry == {"M1_e": -1, "M1_c": 1}
    assert "r1" in map(lambda p: p.enzyme_id, mi.priors.kcat_priors)
    assert all(c in experiment_ids for c in expected_experiment_ids)
    for key, expected_codes in expected_stan_codes.items():
        assert mi.stan_codes.__dict__[key] == expected_codes
