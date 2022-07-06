"""Unit tests for model ode function."""

import sys


if sys.version_info[1] <= 9:
    from importlib_resources import files
else:
    from importlib.resources import files

import cmdstanpy
import pytest
from numpy import isclose


DATA_PATH = files("maud.data")
MODEL_PATH = files("maud.stan").joinpath("model.stan")

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
        [
            "conc[1,1]",
            "conc[1,2]",
            "conc[1,3]",
            "conc[1,4]",
            "flux[1,1]",
            "flux[1,4]",
        ],
        ["conc|A", "conc|B", "conc|C", "conc|D", "flux|r1", "flux|r4"],
        [5.0, 0.323117, 3.02187, 0.5, 0.421816, 2.11674],
    ),
    (
        "methionine",
        [
            "conc[1,5]",
            "flux[3,1]",
            "allostery[4,5]",
        ],
        [
            "conc|amet_c|dataset1",
            "flux|METAT|dataset6",
            "allostery|GNMT1_GNMT|dataset9",
        ],
        [5.931050e-5, 0.00108749, 0.0623294],
    ),
]


@pytest.mark.parametrize(
    "folder,csv_ixs,meaningful_ixs,expected_values", test_cases
)
def test_model_ode(folder, csv_ixs, meaningful_ixs, expected_values):
    """Test that the function get_input_data behaves as expected."""
    input_path = DATA_PATH / folder
    init_data = input_path / "inits.json"
    input_data = input_path / "input_data_train.json"
    SIM_CONFIG["inits"] = init_data
    model = cmdstanpy.CmdStanModel(stan_file=MODEL_PATH)
    sim_values = (
        model.sample(data=input_data.__str__(), **SIM_CONFIG)
        .draws_pd()
        .loc[0, csv_ixs]
        .to_list()
    )
    for meaningful_ix, expected_value, sim_value in zip(
        meaningful_ixs, expected_values, sim_values
    ):
        var = meaningful_ix.split("|")[0]
        ix = " ".join(meaningful_ix.split("|")[1:])
        msg = (
            f"\nExpected value of {var} {ix}: {expected_value}"
            f"\nSimulated value:\n {sim_value}"
        )
        assert isclose(expected_value, sim_value), msg
