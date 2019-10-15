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
