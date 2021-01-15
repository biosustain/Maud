"""Script for running a simulation study."""

import json
import os
from copy import deepcopy

import numpy as np
from cmdstanpy import CmdStanMCMC, CmdStanModel
from cmdstanpy.utils import jsondump
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
    "abs_tol": 1e-7,
    "rel_tol": 1e-7,
    "max_num_steps": int(1e9),
    "likelihood": 1,
    "timepoint": 500,
}
SIM_CONFIG = {
    "chains": 1,
    "iter_sampling": 1,
    "fixed_param": True,
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


def add_measurements_to_maud_input(
    mi: MaudInput, sim: CmdStanMCMC, input_data: dict
) -> MaudInput:
    """Replace the measurements in a maud input with simulated ones."""
    out = deepcopy(mi)
    code_to_exp = {v: k for k, v in mi.stan_codes.experiment_codes.items()}
    code_to_enz = {v: k for k, v in mi.stan_codes.enzyme_codes.items()}
    code_to_mic = {v: k for k, v in mi.stan_codes.mic_codes.items()}
    code_to_rxn = {v: k for k, v in mi.stan_codes.reaction_codes.items()}
    var_ids = {
        "yenz_sim": [code_to_enz[i] for i in input_data["enzyme_yenz"]],
        "yconc_sim": [code_to_mic[i] for i in input_data["mic_ix_yconc"]],
        "yflux_sim": [code_to_rxn[i] for i in input_data["reaction_yflux"]],
    }
    exp_ids = {
        "yenz_sim": [code_to_exp[i] for i in input_data["experiment_yenz"]],
        "yconc_sim": [code_to_exp[i] for i in input_data["experiment_yconc"]],
        "yflux_sim": [code_to_exp[i] for i in input_data["experiment_yflux"]],
    }
    for var, measurement_type in zip(
        ["yenz_sim", "yconc_sim", "yflux_sim"], ["enzyme", "metabolite", "reaction"]
    ):
        if var in sim.stan_variable_dims.keys():
            for i, y in enumerate(sim.stan_variable(var).values[0]):
                exp_id, var_id = exp_ids[var][i], var_ids[var][i]
                experiment = [e for e in out.experiments.experiments if e.id == exp_id][
                    0
                ]
                measurement = experiment.measurements[measurement_type][var_id]
                measurement.value = y
    return out


def enrich_true_values(tvin, input_data):
    """Add true values for auxiliary parameters."""

    def logz_for_vec(truth, loc, scale):
        t = np.array(truth)
        lc = np.array(loc)
        s = np.array(scale)
        return (np.log(t) - np.log(lc)) / s

    def z_for_vec(truth, loc, scale):
        t = np.array(truth)
        lc = np.array(loc)
        s = np.array(scale)
        return (t - lc) / s

    def logz_for_mat(truth, loc, scale):
        return np.array(
            [
                (np.log(np.array(loc[i])) - np.log(np.array(truth[i])))
                / np.array(scale[i])
                for i, _ in enumerate(truth)
            ]
        )

    def z_for_mat(truth, loc, scale):
        return np.array(
            [
                (np.array(loc[i]) - np.array(truth[i])) / np.array(scale[i])
                for i, _ in enumerate(truth)
            ]
        )

    return {
        **tvin,
        **{
            "drain_z": z_for_vec(
                tvin["drain"],
                input_data["prior_loc_drain"],
                input_data["prior_scale_drain"],
            ),
            "log_km_z": logz_for_vec(
                tvin["km"], input_data["prior_loc_km"], input_data["prior_scale_km"]
            ),
            "log_kcat_z": logz_for_vec(
                tvin["kcat"],
                input_data["prior_loc_kcat"],
                input_data["prior_scale_kcat"],
            ),
            "log_ki_z": logz_for_vec(
                tvin["ki"], input_data["prior_loc_ki"], input_data["prior_scale_ki"]
            ),
            "log_dissociation_constant_t_z": logz_for_vec(
                tvin["dissociation_constant_t"],
                input_data["prior_loc_diss_t"],
                input_data["prior_scale_diss_t"],
            ),
            "log_dissociation_constant_r_z": logz_for_vec(
                tvin["dissociation_constant_r"],
                input_data["prior_loc_diss_r"],
                input_data["prior_scale_diss_r"],
            ),
            "log_transfer_constant_z": logz_for_vec(
                tvin["transfer_constant"],
                input_data["prior_loc_tc"],
                input_data["prior_scale_tc"],
            ),
            "log_enzyme_z": logz_for_mat(
                tvin["conc_unbalanced"],
                input_data["prior_loc_unbalanced"],
                input_data["prior_scale_unbalanced"],
            ),
            "log_conc_unbalanced_z": logz_for_mat(
                tvin["conc_unbalanced"],
                input_data["prior_loc_unbalanced"],
                input_data["prior_scale_unbalanced"],
            ),
            "log_drain_z": z_for_mat(
                tvin["conc_unbalanced"],
                input_data["prior_loc_unbalanced"],
                input_data["prior_scale_unbalanced"],
            ),
            "log_phos_kcat_z": logz_for_vec(
                tvin["phos_kcat"],
                input_data["prior_loc_phos_kcat"],
                input_data["prior_scale_phos_kcat"],
            ),
            "log_phos_kcat_z": logz_for_vec(
                tvin["phos_enzyme_conc"],
                input_data["prior_loc_phos_conc"],
                input_data["prior_scale_phos_conc"],
            ),
        },
    }


def main():
    """Run the script."""
    print("Compiling Stan model...")
    model = CmdStanModel(stan_file=STAN_PROGRAM_PATH)

    print("Generating simulation input...")
    mi_sim = load_maud_input_from_toml(TOML_PATH)
    input_data_sim = get_input_data(mi_sim, **ODE_CONFIG)
    with open(TRUE_PARAM_PATH, "r") as f:
        print(TRUE_PARAM_PATH)
        true_values_in = json.load(f)
    true_values = enrich_true_values(true_values_in, input_data_sim)

    print("Simulating...")
    sim = model.sample(data=input_data_sim, inits=true_values, **SIM_CONFIG)

    print("Generating fit input...")
    mi = add_measurements_to_maud_input(mi_sim, sim, input_data_sim)
    input_data_fit = get_input_data(mi, **ODE_CONFIG)
    jsondump("input_data.json", input_data_fit)

    print("Fitting...")
    model.sample(data=input_data_fit, inits=true_values, **FIT_CONFIG)

    print("Analysing results...")
    with open(TRUE_PARAM_PATH, "r") as f:
        true_params = json.load(f)
    csvs = [os.path.join(HERE, f) for f in os.listdir(HERE) if f.endswith(".csv")]
    infd = load_infd(csvs, mi)
    exp_codes = mi.stan_codes.experiment_codes
    enz_codes = mi.stan_codes.enzyme_codes
    met_codes = mi.stan_codes.metabolite_codes
    mic_codes = mi.stan_codes.mic_codes
    unbalanced_mics = input_data_fit["unbalanced_mic_ix"]
    unb_codes = {k: v for k, v in mic_codes.items() if v in unbalanced_mics}
    km_ids = infd.posterior.coords["km_id"].to_series().tolist()
    km_codes = dict(zip(km_ids, range(1, len(km_ids) + 1)))
    true_enz = {
        exp_name: dict(zip(enz_codes.keys(), true_params["enzyme"][i]))
        for i, exp_name in enumerate(exp_codes.keys())
    }
    true_unb = {
        exp_name: dict(zip(unb_codes.keys(), true_params["conc_unbalanced"][i]))
        for i, exp_name in enumerate(exp_codes.keys())
    }
    yconc = input_data_fit["yconc"]
    exp_yconc = [
        list(exp_codes.keys())[i - 1] for i in input_data_fit["experiment_yconc"]
    ]
    mic_yconc = [list(mic_codes.keys())[i - 1] for i in input_data_fit["mic_ix_yconc"]]
    meas_conc = {
        exp_name: {
            mic_yconc[i]: yconc[i] for i, e in enumerate(exp_yconc) if e == exp_name
        }
        for exp_name in exp_yconc
    }
    for varname, codes, figsize, truth, meas in (
        ["enzyme", enz_codes, [16, 6], true_enz, None],
        ["conc", unb_codes, [16, 6], true_unb, meas_conc],
    ):
        f, axes = plot_experiment_var(
            infd, varname, codes, exp_codes, truth, meas, logscale=True
        )
        f.set_size_inches(figsize)
        f.savefig(f"{varname}_posteriors.png", bbox_inches="tight")
    for varname, codes, figsize, logscale in zip(
        ["formation_energy", "km", "kcat"],
        [met_codes, km_codes, enz_codes],
        [[15, 8], [15, 10], [15, 8]],
        [False, True, True],
    ):
        f, _ = plot_1d_var(
            infd, varname, codes, true_params[varname], logscale=logscale
        )
        plt.tight_layout()
        f.set_size_inches(figsize)
        f.savefig(f"{varname}_posteriors.png", bbox_inches="tight")
    plt.close("all")

    print(f"Done! See {OUTPUT_DIR_FIT} for output csvs and plots.")


if __name__ == "__main__":
    main()
