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
TOML_PATH = os.path.join(HERE, "../tests/data/ecoli_small_experiments.toml")
TRUE_PARAMS = {
    "km": [
        3.0, 0.16, 0.0125, 0.06, 15.0, 0.025, 0.03, 0.02, 15.0, 0.025, 1.7,
        0.69, 1.0, 0.5, 0.88, 0.1, 0.451, 2.4, 2.0, 0.33, 0.079
    ],
    "ki": [0.26],
    "dissociation_constant_t": [1, 1],
    "dissociation_constant_r": [0.01],
    "transfer_constant": [1, 1],
    "kcat": [1550.0, 110.0, 56.0, 24.0, 5.67, 4.13, 6.0],
    "conc_unbalanced": [
        [2.080410885611594, 0.6113649750632913, 5.408003186466161, 0.1, 0.1, 0.24187807672413633],
        [14.302230807530746, 0.9898468738819469, 5.669676161196416, 0.1, 0.1, 0.33726831321487066],
        [6.354952226318044, 0.7113810622662209, 2.843811705692167, 0.1, 0.1, 1.692546934887993],
        [1.8990486465075278, 0.9105422318194938, 3.8454204039099675, 0.1, 0.1, 0.2833704100025282],
        [2.58774462726839, 0.9806673979313613, 3.8811478289600676, 0.1, 0.1, 0.20534440886003383]
    ],
    "enzyme": [
        [0.03338748587758992, 0.03233687344889826, 0.00470561182567603, 0.0057128462581434464, 0.0704592675242211, 0.004038784112975896, 0.019814776219274497],
        [0.1, 0.013087078730074364, 0.00365161452395729, 0.01173326160905437, 0.05228843722116154, 0.008492123537417372, 0.022968874089049705],
        [0.044449648932757366, 0.04634699167660036, 0.008006922986439756, 0.04362849543945923, 0.06197117681867592, 0.028518936258079924, 0.020866940984252118],
        [0.03421127030955773, 0.04930313430995316, 0.007853203230449156, 0.0075998055009798265, 0.06783806498863625, 0.0025627145433086006, 0.02367223995109299],
        [0.053645551520857676, 0.04484149627257232, 0.003917011890185323, 0.006854744097845603, 0.06529811635401973, 0.009351216672916836, 0.1]
    ],
    "formation_energy": [
        -1336.3, -1333.8, -2220.9, -1440.8, -2313.0, -1073.3, -1111.9,
        -1106.4, 0.0
    ],
    "drain": [[-1.539], [0.563], [-0.456], [0.194], [0.1376979321]]
}
ODE_CONFIG = {
    "abs_tol": 1e-6,
    "rel_tol": 1e-6,
    "max_num_steps": int(1e9),
    "likelihood": 1,
    "timepoint": 500,
}
SIM_CONFIG = {
    "chains": 1,
    "iter_sampling": 1,
    "output_dir": OUTPUT_DIR_SIM,
    "fixed_param": True
}
FIT_CONFIG = {
    "chains": 4,
    "parallel_chains": 2,
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


def main():
    print("Compiling Stan model...")
    model = CmdStanModel(stan_file=STAN_PROGRAM_PATH)

    print("Generating simulation input...")
    mi_sim = load_maud_input_from_toml(TOML_PATH)
    input_data_sim = get_input_data(mi_sim, **ODE_CONFIG)
    true_params = TRUE_PARAMS.copy()
    true_params["formation_energy_z"] = standardise_list(
        true_params["formation_energy"],
        input_data_sim["prior_loc_formation_energy"],
        input_data_sim["prior_scale_formation_energy"]
    )

    print("Simulating...")
    sim = model.sample(data=input_data_sim, inits=true_params, **SIM_CONFIG)

    print("Generating fit input...")
    mi_fit = add_measurements_to_maud_input(mi_sim, sim, input_data_sim)
    input_data_fit = get_input_data(mi_fit, **ODE_CONFIG)
    # inits_fit = get_initial_conditions(input_data_fit, mi_fit)
    inits_fit = {k: list(v) * FIT_CONFIG["chains"] for k, v in TRUE_PARAMS.items()}

    print("Fitting...")
    fit = model.sample(data=input_data_fit, inits=true_params, **FIT_CONFIG)

    print(f"Done! See {output_dir} for output csvs.")

    
if __name__ == "__main__":
    main()
