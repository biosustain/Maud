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


def test_create_haldane_line():
    """Test that the function create_haldane_line behaves as expected."""
    rxn_id = "tr"
    mechanism = "ordered_unibi"
    param_codes = {"tr_Keq": 1, "tr_Kia": 2, "tr_Kq": 3, "tr_Kcat1": 4, "tr_Kcat2": 5}
    expected_line = "real tr_Kip = get_Kip_ordered_unibi(p[1],p[2],p[3],p[4],p[5]);"
    generated_line = code_generation.create_haldane_line(
        param_codes, "Kip", mechanism, rxn_id, code_generation.HALDANE_PARAMS
    )
    assert generated_line == expected_line


def test_get_regulatory_string():
    """Test that the function get_regulatory_string behaves as expected."""
    inhibitor_codes = {"m1": 1, "m2": 2}
    param_codes = {
        "e1_dissociation_constant_t_m1": 2,
        "e1_dissociation_constant_t_m2": 3,
        "e1_transfer_constant": 4,
    }
    enzyme_name = "e1"
    generated = code_generation.get_regulatory_string(
        inhibitor_codes, param_codes, enzyme_name
    )
    expected = (
        "get_regulatory_effect("
        "empty_array,{m[1],m[2]},free_enzyme_ratio_e1,empty_array,{p[2],p[3]},p[4]"
        ")"
    )
    assert generated == expected


def test_mechanism_templates():
    """Test that mechanism templates are as expected"""
    codes = {
        "S0": 1,
        "S1": 2,
        "S2": 3,
        "P0": 4,
        "P1": 5,
        "enz": 1,
        "Keq": 2,
        "Kcat1": 3,
        "Kcat2": 4,
        "Ka": 5,
        "Kb": 6,
        "Kc": 7,
        "Kp": 8,
        "Kq": 9,
        "Kia": 10,
        "Kib": 11,
        "Kib": 12,
        "Kiq": 13,
        "Tr": 1,
        "Dr": 2,
    }
    enzyme_id = "e1"
    expected_calls = {
        "uniuni": "uniuni(m[1],m[4],p[1]*p[3],p[1]*p[4],p[5],p[2])",
        "ordered_unibi": (
            "ordered_unibi(m[1],m[4],m[5],p[1]*p[3],p[1]*p[4],p[5],"
            "p[8],p[9],p[10],e1_Kip,e1_Kiq],p[2])"
        ),
        "ordered_bibi": (
            "ordered_bibi(m[1],m[2],m[4],m[5],p[1]*p[3],p[1]*p[4],p[5],"
            "p[6],p[8],p[9],e1_Kia,p[12],e1_Kip,p[13],p[2])"
        ),
        "pingpong": (
            "pingpong(m[1],m[2],m[4],m[5],p[1]*p[3],p[1]*p[4],p[5],p[6],"
            "p[8],p[9],p[10],p[12],e1_Kip,p[13],p[2])"
        ),
        "ordered_terbi": (
            "ordered_terbi(m[1],m[2],m[3],m[4],m[5],p[1]*p[3],p[1]*p[4],"
            "p[5],p[6],p[7],e1_Kp,p[9],p[10],p[12],p[],e1_Kip,p[13],p[2])"
        ),
        "modular_rate_law": "modular_rate_law(1,2)",
    }
    generated_calls = {
        mechanism: template.render({**codes, **{"enz_id": enzyme_id}})
        for mechanism, template in code_generation.MECHANISM_TEMPLATES.items()
    }
    for mechanism, generated_call in generated_calls.items():
        expected_call = expected_calls[mechanism]
        assert generated_call == expected_call
