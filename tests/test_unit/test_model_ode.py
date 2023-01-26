"""Unit tests for model ode function."""

import os

import cmdstanpy
import pytest
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
    "output_dir": None,
}

test_cases = [
    (
        "example_ode",
        {
            "conc_train": {
                "B": ((0, 0), 5.0),
                "C": ((0, 1), 0.323117),
                "D": ((0, 2), 3.02187),
            },
            "flux_train": {
                "r1": ((0, 0), 0.421816),
                "r4": ((0, 3), 2.11674),
            },
        },
    ),
    #  methionine input commented out until we find the right values
    # (
    #     "methionine",
    #     {
    #         "conc_train": {
    #             "amet_c|dataset1": ((1, 5), 5.931050e-5)
    #         },
    #         "flux_train": {
    #             "METAT|dataset6": ((3, 1), 0.00108749),
    #         },
    #         "allostery": {
    #             "GNMT1_GNMT|dataset9": ((4, 5), 0.0623294)
    #         }
    #     },
    # ),
]


@pytest.mark.parametrize("folder,expected_values", test_cases)
def test_model_ode(folder, expected_values):
    """Test that the function get_input_data behaves as expected."""
    input_path = os.path.join(data_path, folder)
    init_data = os.path.join(input_path, "inits.json")
    input_data = os.path.join(input_path, "input_data_train.json")
    SIM_CONFIG["inits"] = init_data
    model = cmdstanpy.CmdStanModel(stan_file=model_path)
    mcmc = model.sample(data=input_data, **SIM_CONFIG)
    msg = (
        "\nExpected value of {var} {ix}: {expected_value}"
        "\nSimulated value:\n {sim_value}"
    )
    for var, expected in expected_values.items():
        sim_values = mcmc.draws_xr(vars=var).sel(chain=1, draw=0)[var]
        for _, (ix, expected_value) in expected.items():
            sim_value = sim_values.values[ix]
            assert isclose(expected_value, sim_value), msg.format(
                var=var, ix=ix, expected_value=expected_value
            )
