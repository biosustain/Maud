"""Script for running a simulation study based on the data at """

from maud.io import load_maud_input_from_toml
from maud.sampling import get_input_data, get_initial_conditions
from cmdstanpy import CmdStanModel, CmdStanMCMC
from cmdstanpy.utils import jsondump
from maud.data_model import MaudInput
from copy import deepcopy
import numpy as np
import pandas as pd
import os

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
    "output_dir": OUTPUT_DIR_SIM,
    "fixed_param": True
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


def add_measurements_to_maud_input(mi: MaudInput, sim: CmdStanMCMC, input_data: dict):
    out = deepcopy(mi)
    code_to_exp = {v: k for k, v in mi.stan_codes.experiment_codes.items()}
    code_to_enz = {v: k for k, v in mi.stan_codes.enzyme_codes.items()}
    code_to_mic = {v: k for k, v in mi.stan_codes.mic_codes.items()}
    code_to_rxn = {v: k for k, v in mi.stan_codes.reaction_codes.items()}
    yenz_experiments = [code_to_exp[i] for i in  input_data["experiment_yenz"]]
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


def standardise_list(l, mean, scale):
    return ((np.array(l) - np.array(mean)) / np.array(scale)).tolist()


def get_2d_var_df(mcmc, varname, ix, cols):
    vals = mcmc.stan_variable(varname).values
    dims = mcmc.stan_variable_dims[varname]
    return pd.DataFrame(vals.reshape(dims, order="F"), index=ix, columns=cols)


def main():
    print("Compiling Stan model...")
    model = CmdStanModel(stan_file=STAN_PROGRAM_PATH)

    print("Generating simulation input...")
    mi_sim = load_maud_input_from_toml(TOML_PATH)
    input_data_sim = get_input_data(mi_sim, **ODE_CONFIG)

    print("Simulating...")
    sim = model.sample(data=input_data_sim, **SIM_CONFIG)

    print("Generating fit input...")
    mi_fit = add_measurements_to_maud_input(mi_sim, sim, input_data_sim)
    input_data_fit = get_input_data(mi_fit, **ODE_CONFIG)

    print("Fitting...")
    fit = model.sample(data=input_data_sim, **FIT_CONFIG)

    print(f"Done! See {output_dir} for output csvs.")

    
if __name__ == "__main__":
    main()
