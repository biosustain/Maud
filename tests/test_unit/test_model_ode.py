"""Unit tests for model ode function."""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

import importlib_resources
import pytest
from numpy import isclose

from maud.data.example_inputs import example_ode
from maud.running_stan import load_stan_model


@dataclass
class ExpectedValue:
    """An expected value."""

    name: str
    ix: Tuple[int, int]
    value: float


@dataclass
class ODETestCase:
    """A test case."""

    name: str
    inits_path: Path
    input_data_path: Path
    expected_conc_train: List[ExpectedValue]
    expected_flux_train: List[ExpectedValue]


ODE_CONFIG_OPTIONS = [
    "input_data_train_adjoint_algebra.json",
    "input_data_train_adjoint.json",
    "input_data_train_bdf_algebra.json",
    "input_data_train_bdf.json",
]

TEST_CASES = [
    ODETestCase(
        name="example_ode",
        inits_path=importlib_resources.files(example_ode).joinpath(
            "inits.json"
        ),
        input_data_path=importlib_resources.files(example_ode).joinpath(
            varied_ode_config
        ),
        expected_conc_train=[
            ExpectedValue("B", (0, 0), 5.0),
            ExpectedValue("C", (0, 1), 0.323117),
            ExpectedValue("D", (0, 2), 3.02187),
        ],
        expected_flux_train=[
            ExpectedValue("r1", (0, 0), 0.421816),
            ExpectedValue("r4", (0, 3), 2.11674),
        ],
    )
    for varied_ode_config in ODE_CONFIG_OPTIONS
]


@pytest.mark.parametrize("test_case", TEST_CASES)
def test_model_ode(test_case):
    """Test that the function get_input_data behaves as expected."""
    model = load_stan_model("model", stanc_options={}, cpp_options={})
    mcmc = model.sample(
        data=str(test_case.input_data_path),
        inits=str(test_case.inits_path),
        chains=1,
        fixed_param=True,
        iter_warmup=0,
        iter_sampling=1,
        show_progress=False,
        show_console=True,
        threads_per_chain=1,
        output_dir=None,
    )
    msg = (
        "\nExpected value of {var} {ix}: {expected_value}"
        "\nSimulated value:\n {sim_value}"
    )
    for var in ["conc_train", "flux_train"]:
        sim_values = mcmc.draws_xr(vars=[var]).sel(chain=1, draw=0)[var].values
        for expected_value in getattr(test_case, "expected_" + var):
            sim_value = sim_values[expected_value.ix]
            assert isclose(expected_value.value, sim_value), msg.format(
                var=var,
                ix=expected_value.ix,
                expected_value=expected_value.value,
                sim_value=sim_value,
            )
