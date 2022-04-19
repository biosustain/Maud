"""Unit tests for model ode function."""

import os

import cmdstanpy
from numpy import isclose


here = os.path.dirname(__file__)
data_path = os.path.join(here, "..", "data")
model_path = os.path.join(here, "..", "..", "src", "maud", "stan", "model.stan")

SIM_CONFIG = {
    "chains": 1,
    "fixed_param": True,
    "iter_warmup": 0,
    "iter_sampling": 1,
    "show_progress": False,
    "threads_per_chain": 1,
}


def test_model_ode():
    """Test that the function get_input_data behaves as expected."""
    input_path = os.path.join(data_path, "example_ode")
    init_data = os.path.join(input_path, "inits.json")
    input_data = os.path.join(input_path, "input_data.json")
    SIM_CONFIG["inits"] = init_data
    model = cmdstanpy.CmdStanModel(stan_file=model_path)
    remap = {
        "conc[1,1]": "A",
        "conc[1,2]": "B",
        "conc[1,3]": "C",
        "conc[1,4]": "D",
        "flux[1,1]": "r1",
        "flux[1,4]": "r4",
    }
    true_values = {
        "A": 5.0,
        "B": 0.323117,
        "C": 3.02187,
        "D": 0.5,
        "r1": 0.421816,
        "r4": 2.11674,
    }
    sim_values = (
        model.sample(data=input_data, **SIM_CONFIG)
        .draws_pd()
        .rename(columns=remap)
        .T.loc[true_values.keys(), 0]
        .to_dict()
    )
    msg = f"\nTrue values:\n {true_values}\nSimulated values:\n {sim_values}"
    for true_value, sim_value in zip(
        true_values.values(), sim_values.values()
    ):
        assert isclose(true_value, sim_value), msg
