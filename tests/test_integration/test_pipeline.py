"""Integration test from code generation to correct sampling."""

import os
import shutil
import tempfile

import arviz as az
import pandas as pd

import maud.io as io
import maud.sampling as sampling


here = os.path.dirname(__file__)
data_path = os.path.join(here, "../data")


def get_keys(stan_codes):
    """Retrive keys from stan codes."""
    return [key for key in stan_codes.keys()]


def get_measurement_index(mi, mic_codes, enzyme_codes, experiment_codes):
    """Return measurement index to align log-likelihoods."""
    measurement_index = [
        f"{meas_id}//{exp_name}//{value.value}"
        for exp_name, exp in mi.experiments.items()
        for meas_type, measurement in exp.measurements.items()
        for meas_id, value in measurement.items()
        if meas_type not in ["enzyme"]
    ]

    measurement_index_ordered = pd.DataFrame()
    measurement_index_ordered["measurement_id"] = measurement_index
    measurement_index_ordered_split = measurement_index_ordered[
        "measurement_id"
    ].str.split("//")

    measurement_data = pd.DataFrame(
        measurement_index_ordered_split.to_list(),
        columns=["measurement_id", "experiment", "value"],
    )
    measurement_index_new = measurement_index_ordered.copy(deep=True)
    measurement_index_new["measurement_id"] = measurement_data["measurement_id"]
    measurement_index_new["experiment"] = measurement_data["experiment"]
    measurement_index_new["value"] = measurement_data["value"]

    meas_codes = mic_codes + enzyme_codes
    meas_index = pd.DataFrame(meas_codes, columns=["measurement_id"])
    meas_index["order"] = [i for i in range(meas_index.shape[0])]

    measurement_index_new = measurement_index_new.merge(
        meas_index, on="measurement_id", how="left"
    )

    for exp in experiment_codes:
        meas_idx_temp = (
            measurement_index_new[measurement_index_new["experiment"] == exp]
            .sort_values(by="order")
            .copy(deep=True)
            .values
        )
        measurement_index_new[
            measurement_index_new["experiment"] == exp
        ] = meas_idx_temp

    return measurement_index_new


def get_coords_dims(mi):
    """Return coordinates for models."""
    experiment_codes = get_keys(mi.stan_codes["experiment"])
    reaction_codes = get_keys(mi.stan_codes["reaction"])
    enzyme_codes = get_keys(mi.stan_codes["enzyme"])
    metabolite_codes = get_keys(mi.stan_codes["metabolite"])
    kinetic_parameter_codes = get_keys(mi.stan_codes["kinetic_parameter"])
    mic_codes = get_keys(mi.stan_codes["metabolite_in_compartment"])

    measurement_index = get_measurement_index(
        mi=mi,
        mic_codes=mic_codes,
        enzyme_codes=enzyme_codes,
        experiment_codes=experiment_codes,
    )

    coords = {
        "reactions": reaction_codes,
        "mic": mic_codes,
        "metabolites": metabolite_codes,
        "experiments": experiment_codes,
        "parameter_names": kinetic_parameter_codes,
        "measurement_names": measurement_index["measurement_id"]
        + "//"
        + measurement_index["experiment"],
    }

    dims = {
        "conc": ["experiments", "mic"],
        "flux": ["experiments", "reactions"],
        "kinetic_parameters": ["parameter_names"],
        "delta_g": ["reactions"],
        "formation_energy": ["metabolites"],
        "log_like": ["measurement_names"],
    }

    return coords, dims


def test_linear():
    """
    Tests linear model.

    tests from code generation to sampling of the linear model by computing 200
    samples after 200 warmups and checking if the sampled median is within the
    94% CI of a precomputed set.
    """

    linear_input = os.path.join(data_path, "linear.toml")
    temp_directory = tempfile.mkdtemp(dir=data_path)
    control_directory = os.path.join(data_path, "linear_control_set")
    linear_control_files = [
        os.path.join(control_directory, f"inference_model_linear_test-{i}.csv")
        for i in range(1, 5)
    ]

    mi = io.load_maud_input_from_toml(linear_input)

    coords, dims = get_coords_dims(mi)

    linear_input_values = {
        "f_tol": 1e-6,
        "rel_tol": 1e-9,
        "max_steps": int(1e9),
        "likelihood": 1,
        "n_samples": 200,
        "n_warmup": 200,
        "n_chains": 4,
        "n_cores": 4,
        "timepoint": 500,
        "output_dir": temp_directory,
        "data_path": linear_input,
    }

    linear_fit_test = sampling.sample(**linear_input_values)
    linear_fit_control = az.from_cmdstan(
        posterior=linear_control_files, coords=coords, dims=dims
    )

    interval_test = pd.DataFrame()

    interval_test["sample_mean"] = az.summary(linear_fit_test)["mean"].values
    interval_test["lower"] = az.summary(linear_fit_control)["hpd_3%"].values
    interval_test["upper"] = az.summary(linear_fit_control)["hpd_97%"].values
    interval_test["coord"] = az.summary(linear_fit_test).index.values
    interval_test["within"] = float("NaN")
    interval_test.set_index("coord", inplace=True)

    for i in range(0, interval_test.shape[0]):
        if (
            interval_test["lower"][i]
            <= interval_test["sample_mean"][i]
            <= interval_test["upper"][i]
        ):
            interval_test["within"][i] = 1
        else:
            interval_test["within"][i] = 0

    out_of_bounds = interval_test[interval_test["within"] == 0]

    shutil.rmtree(temp_directory)

    assert out_of_bounds.empty, out_of_bounds


test_linear()
