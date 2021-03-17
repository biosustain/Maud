"""Test that simple simulation studies give expected output."""

import json
import os

import numpy as np
import pytest

from maud.analysis import load_infd
from maud.io import load_maud_input_from_toml
from maud.simulation_study import run_simulation_study


HERE = os.path.dirname(__file__)
CASES_DIR = os.path.join(HERE, "..", "data", "example_features")
MAUD_PATH = os.path.join(HERE, "..", "..", "src", "maud")
CASES = map(
    lambda f: os.path.join(CASES_DIR, f),
    [
        "michaelis_menten",
        "drain",
        "phosphorylation",
        "one_allosteric_modifier",
        "two_allosteric_modifiers",
    ],
)
TRUE_PARAMS_FILENAME = "true_params.json"


@pytest.mark.parametrize("input_dirname", CASES)
def test_linear(input_dirname):
    """Test that the linear model works."""

    input_dir_path = os.path.join(HERE, input_dirname)
    mi_in = load_maud_input_from_toml(input_dir_path)
    true_params_path = os.path.join(input_dir_path, TRUE_PARAMS_FILENAME)
    with open(true_params_path, "r") as f:
        true_params = json.load(f)
    study = run_simulation_study(mi_in, true_params)
    infd = load_infd(study.samples.runset.csv_files, study.mi)
    for param_name, param_vals in true_params.items():
        if any(param_vals):
            dimnames = [
                d
                for d in infd.posterior[param_name].dims
                if d not in ["chain", "draw"]
            ]
            q = (
                infd.posterior[param_name]
                .to_series()
                .unstack(dimnames)
                .quantile([0.025, 0.975])
                .T.assign(true=np.array(param_vals).ravel())
            )
            q.columns = ["low", "high", "true"]
            for i, row in q.iterrows():
                msg = (
                    f"True value for {param_name} outside 95% CI at coord {str(i)}!\n"
                    f"\tTrue value: {str(row['true'])}\n"
                    f"\t2.5% posterior quantile: {str(row['low'])}\n"
                    f"\t97.5% posterior quantile: {str(row['high'])}\n"
                )
                assert (
                    row["true"] >= row["low"] and row["true"] <= row["high"]
                ), msg
