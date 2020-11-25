"""Script for running a simulation study."""

import json
import os
from copy import deepcopy

from cmdstanpy import CmdStanMCMC, CmdStanModel
from matplotlib import pyplot as plt

from maud.analysis import load_infd, plot_1d_var, plot_experiment_var
from maud.data_model import MaudInput
from maud.io import load_maud_input_from_toml
from maud.sampling import get_input_data


HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR_SIM = HERE
OUTPUT_DIR_FIT = HERE
STAN_PROGRAM_PATH = os.path.join(HERE, "../src/maud/inference_model.stan")
TOML_PATH = os.path.join(HERE, "../tests/data/ecoli_small.toml")
TRUE_PARAM_PATH = os.path.join(HERE, "../tests/data/true_params_ecoli_small.json")
ODE_CONFIG = {
    "abs_tol": 1e-6,
    "rel_tol": 1e-6,
    "max_num_steps": int(1e9),
    "likelihood": 1,
    "timepoint": 500,
}
SIM_CONFIG = {
    "chains": 1,
    "inits": TRUE_PARAM_PATH,
    "iter_sampling": 1,
    "fixed_param": True,
}
FIT_CONFIG = {
    "chains": 4,
    "parallel_chains": 2,
    "inits": TRUE_PARAM_PATH,
    "iter_warmup": 200,
    "iter_sampling": 200,
    "output_dir": OUTPUT_DIR_FIT,
    "max_treedepth": 11,
    "save_warmup": True,
    "show_progress": True,
    "step_size": 0.025,
    "adapt_delta": 0.99,
}


def add_measurements_to_maud_input(
    mi: MaudInput, sim: CmdStanMCMC, input_data: dict
) -> MaudInput:
    """Replace the measurements in a maud input with simulated ones."""
    out = deepcopy(mi)
    code_to_exp = {v: k for k, v in mi.stan_codes.experiment_codes.items()}
    code_to_enz = {v: k for k, v in mi.stan_codes.enzyme_codes.items()}
    code_to_mic = {v: k for k, v in mi.stan_codes.mic_codes.items()}
    code_to_rxn = {v: k for k, v in mi.stan_codes.reaction_codes.items()}
    yenz_experiments = [code_to_exp[i] for i in input_data["experiment_yenz"]]
    yenz_enzymes = [code_to_enz[i] for i in input_data["enzyme_yenz"]]
    yconc_experiments = [code_to_exp[i] for i in input_data["experiment_yconc"]]
    yconc_mics = [code_to_mic[i] for i in input_data["mic_ix_yconc"]]
    yflux_experiments = [code_to_exp[i] for i in input_data["experiment_yflux"]]
    yflux_reactions = [code_to_rxn[i] for i in input_data["reaction_yflux"]]
    for i, y in enumerate(sim.stan_variable("yenz_sim").values[0]):
        exp_id, enz_id = yenz_experiments[i], yenz_enzymes[i]
        experiment = [e for e in out.experiments.experiments if e.id == exp_id][0]
        measurement = experiment.measurements["enzyme"][enz_id]
        measurement.value = y
    for i, y in enumerate(sim.stan_variable("yconc_sim").values[0]):
        exp_id, mic_id = yconc_experiments[i], yconc_mics[i]
        experiment = [e for e in out.experiments.experiments if e.id == exp_id][0]
        measurement = experiment.measurements["metabolite"][mic_id]
        measurement.value = y
    for i, y in enumerate(sim.stan_variable("yflux_sim").values[0]):
        exp_id, rxn_id = yflux_experiments[i], yflux_reactions[i]
        experiment = [e for e in out.experiments.experiments if e.id == exp_id][0]
        measurement = experiment.measurements["reaction"][rxn_id]
        measurement.value = y
    return out


def main():
    """Run the script."""
    print("Compiling Stan model...")
    model = CmdStanModel(stan_file=STAN_PROGRAM_PATH)

    print("Generating simulation input...")
    mi_sim = load_maud_input_from_toml(TOML_PATH)
    input_data_sim = get_input_data(mi_sim, **ODE_CONFIG)

    print("Simulating...")
    sim = model.sample(data=input_data_sim, **SIM_CONFIG)

    print("Generating fit input...")
    mi = add_measurements_to_maud_input(mi_sim, sim, input_data_sim)
    input_data_fit = get_input_data(mi, **ODE_CONFIG)

    print("Fitting...")
    fit = model.sample(data=input_data_sim, **FIT_CONFIG)

    print("Analysing results...")
    with open(TRUE_PARAM_PATH, "r") as f:
        true_params = json.load(f)
    csvs = [f for f in os.listdir(".") if f.endswith(".csv")]
    infd = load_infd(csvs, mi)
    exp_codes = mi.stan_codes.experiment_codes
    enz_codes = mi.stan_codes.enzyme_codes
    met_codes = mi.stan_codes.metabolite_codes
    mic_codes = mi.stan_codes.mic_codes
    unbalanced_mics = input_data_fit["unbalanced_mic_ix"]
    unb_codes = {k: v for k, v in mic_codes.items() if v in unbalanced_mics}
    km_ids = infd.posterior.coords["km_id"].to_series().tolist()
    km_codes = dict(zip(km_ids, range(1, len(km_ids) + 1)))
    for varname, codes, figsize, truth in (
        ["enzyme", enz_codes, [15, 8], true_params["enzyme"]],
        ["conc", unb_codes, [15, 8], true_params["conc_unbalanced"]],
    ):
        f, _ = plot_experiment_var(infd, varname, codes, exp_codes, truth)
        plt.tight_layout()
        f.set_size_inches(figsize)
        f.savefig(f"{varname}_posteriors.png", bbox_inches="tight")
    for varname, codes, figsize in zip(
        ["formation_energy", "km", "kcat"],
        [met_codes, km_codes, enz_codes],
        [[15, 8], [15, 10], [15, 8]],
    ):
        f, _ = plot_1d_var(infd, varname, codes, true_params[varname])
        plt.tight_layout()
        f.set_size_inches(figsize)
        f.savefig(f"{varname}_posteriors.png", bbox_inches="tight")
    plt.close("all")

    print(f"Done! See {OUTPUT_DIR_FIT} for output csvs and plots.")


if __name__ == "__main__":
    main()
