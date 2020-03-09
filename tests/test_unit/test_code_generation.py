"""Unit tests for code generation functions."""

import os

import maud.code_generation as code_generation
import maud.io as io


here = os.path.dirname(__file__)
data_path = os.path.join(here, "../data")


def test_create_stan_program():
    """Test that the function create_stan_program behaves as expected."""
    correct_stan_code_filenames = {"inference": "linear.stan"}
    toml_input_path = os.path.join(data_path, "linear.toml")
    mi = io.load_maud_input_from_toml(toml_input_path)
    for model_type in ["inference"]:
        stan_code_filename = correct_stan_code_filenames[model_type]
        correct_stan_code_path = os.path.join(data_path, stan_code_filename)
        correct_stan_code = open(correct_stan_code_path, "r").read()
        generated_stan_code = code_generation.create_stan_program(mi, model_type)
        correct_lines = correct_stan_code.strip().splitlines(1)
        generated_lines = generated_stan_code.strip().splitlines(1)

        for i, c in enumerate(correct_lines):
            assert c == generated_lines[i]


def test_get_regulatory_string():
    """Test that the function get_regulatory_string behaves as expected."""
    inhibitor_codes = {"m1": "m[1]", "m2": "m[2]"}
    activator_codes = {}
    param_codes = {
        "e1_dissociation_constant_t_m1": 2,
        "e1_dissociation_constant_t_m2": 3,
        "e1_transfer_constant": 4,
    }
    enzyme_name = "e1"
    generated = code_generation.get_regulatory_string(
        inhibitor_codes, activator_codes, param_codes, enzyme_name
    )
    expected = (
        "get_regulatory_effect("
        "empty_array,{m[1],m[2]},free_enzyme_ratio_e1,empty_array,{p[2],p[3]},p[4]"
        ")"
    )
    assert generated == expected


def test_get_modular_rate_codes():
    """Check that the function get_modular_rate_codes works as expected."""
    enz_id = "r"
    competitor_info = []
    substrate_info = [["m1", -1], ["m2", -1]]
    product_info = [["m3", 1], ["m4", 1]]
    par_codes = {"r_keq": 0.1, "r_Ka": 0.1, "r_Kb": 0.3, "r_Kp": 0.2, "r_Kq": 0.2}
    met_codes = {"m1": "m[1]", "m2": "m[2]", "m3": "m[3]", "m4": "m[4]"}
    expected_output = [
        [["m[1]", 0.1, -1], ["m[2]", 0.3, -1]],
        [["m[3]", 0.2, 1], ["m[4]", 0.2, 1]],
        [],
    ]
    assert (
        code_generation.get_modular_rate_codes(
            enz_id, competitor_info, substrate_info, product_info, par_codes, met_codes
        )
        == expected_output
    )
