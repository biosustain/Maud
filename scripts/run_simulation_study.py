"""Run a simulation study using the data at INPUT_DATA."""

import json
import os
import shutil
from datetime import datetime

import arviz as az
import pandas as pd

from maud.cli import sample, simulate
from maud.io import load_maud_input_from_toml


INPUT_DATA = "../../Models/data/march-2021/methionine_cycle"


def get_experiment_table_from_sim(sim_dir: str) -> pd.DataFrame:
    """Get a table of simulated measurements.

    The output should be compatible with maud, so that it is possible to
    overwrite the experiments file.

    """
    ui_path = os.path.join(sim_dir, "user_input")
    sim_csv_path = os.path.join(sim_dir, "samples")
    csv_file = os.path.join(
        sim_csv_path,
        next(filter(lambda f: f.endswith(".csv"), os.listdir(sim_csv_path))),
    )
    with open(os.path.join(sim_csv_path, "input_data.json"), "r") as f:
        stan_input = json.load(f)
    mi = load_maud_input_from_toml(ui_path)
    infd = az.from_cmdstan(csv_file)
    code_to_exp = {v: k for k, v in mi.stan_codes.experiment_codes.items()}
    code_to_mic = {v: k for k, v in mi.stan_codes.mic_codes.items()}
    code_to_rxn = {v: k for k, v in mi.stan_codes.reaction_codes.items()}
    conc_sim = pd.DataFrame(
        {
            "measurement_type": "mic",
            "target_id": map(code_to_mic.get, stan_input["mic_ix_yconc"]),
            "experiment_id": map(code_to_exp.get, stan_input["experiment_yconc"]),
            "measurement": infd.posterior["yconc_sim"].to_series().values,
            "error_scale": stan_input["sigma_conc"],
        }
    )
    flux_sim = pd.DataFrame(
        {
            "measurement_type": "flux",
            "target_id": map(code_to_rxn.get, stan_input["reaction_yflux"]),
            "experiment_id": map(code_to_exp.get, stan_input["experiment_yflux"]),
            "measurement": infd.posterior["yflux_sim"].to_series().values,
            "error_scale": stan_input["sigma_flux"],
        }
    )
    enz_og = pd.read_csv(os.path.join(ui_path, mi.config.experiments_file)).loc[
        lambda df: df["measurement_type"] == "enz"
    ]
    return pd.concat([conc_sim, flux_sim, enz_og], ignore_index=True)


def main():
    """Run the script."""
    here = os.path.dirname(os.path.abspath(__file__))
    now = datetime.now().strftime("%Y%m%d%H%M%S")
    input_path = os.path.join(here, INPUT_DATA)
    input_dirname = os.path.split(input_path)[-1]
    sim_study_folder = os.path.join(here, f"sim_study-{input_dirname}-{now}")
    # make the output directory
    os.mkdir(sim_study_folder)
    # copy the input directory to output/sim_input and output/sample_input
    sim_input_dir, sample_input_dir = (
        shutil.copytree(input_path, os.path.join(sim_study_folder, dirname))
        for dirname in ["sim_input", "sample_input"]
    )
    # simulate some measurements
    sim_dir = simulate(sim_input_dir, output_dir=sim_study_folder, n=1)
    # overwrite the measurements in the sample input based on the simulation
    new_experiments = get_experiment_table_from_sim(sim_dir)
    csv_target = os.path.join(sample_input_dir, load_maud_input_from_toml(sample_input_dir).config.experiments_file)
    new_experiments.to_csv(csv_target)
    # run maud sample against the doctored input
    sample(sample_input_dir, output_dir=sim_study_folder)


if __name__ == "__main__":
    main()
