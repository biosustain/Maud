"""Integration test from code generation to correct sampling."""

import os
import shutil
import tempfile

import pandas as pd

import maud.sampling as sampling


here = os.path.dirname(__file__)
data_path = os.path.join(here, "../data")


def test_linear():
    """
    Tests linear model.

    tests from code generation to sampling of the linear model by computing 200
    samples after 200 warmups and checking if the sampled median is within the
    94% CI of a precomputed set.
    """

    linear_input = os.path.join(data_path, "linear.toml")
    temp_directory = tempfile.mkdtemp(dir=data_path)
    linear_control_file = os.path.join(data_path, "linear_control_samples.csv")
    linear_input_values = {
        "abs_tol": 1e-6,
        "rel_tol": 1e-6,
        "max_num_steps": int(1e9),
        "likelihood": 1,
        "n_samples": 200,
        "n_warmup": 200,
        "n_chains": 1,
        "n_cores": 1,
        "timepoint": 500,
        "output_dir": temp_directory,
        "data_path": linear_input,
        "threads_per_chain": 1,
    }
    fit = sampling.sample(**linear_input_values)
    samples_test = fit.get_drawset()
    samples_ctrl = pd.read_csv(linear_control_file, comment="#").iloc[200:]
    # Check that the output and the control have the same column names.
    assert all(samples_test.columns == samples_ctrl.columns)
    # Check that each output column (other than the diagnostic ones) is
    # statistically similar to its matching control column.
    test_mean = samples_test.mean()
    ctrl_low = samples_ctrl.quantile(0.03)
    ctrl_high = samples_ctrl.quantile(0.97)
    cols = [c for c in samples_test.columns if not c.endswith("__")] + ["lp__"]
    for col in cols:
        assert test_mean[col] >= ctrl_low[col], col + " is too low."
        assert test_mean[col] <= ctrl_high[col], col + " is too high."
    # Delete temporary directory
    shutil.rmtree(temp_directory)
