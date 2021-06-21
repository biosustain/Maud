"""Code for running simulation studies."""


import os
from copy import deepcopy

import numpy as np
from cmdstanpy import CmdStanMCMC, CmdStanModel

from maud.data_model import MaudInput, SimulationStudyOutput
from maud.sampling import STAN_PROGRAM_RELATIVE_PATH, get_input_data


SIM_CONFIG = {"chains": 1, "fixed_param": True}


def add_measurements_to_maud_input(
    mi: MaudInput, sim: CmdStanMCMC, input_data: dict
) -> MaudInput:
    """Replace the measurements in a maud input with simulated ones."""
    out = deepcopy(mi)
    code_to_exp = {i + 1: c for i, c in enumerate(mi.stan_coords.experiments)}
    code_to_enz = {i + 1: c for i, c in enumerate(mi.stan_coords.enzymes)}
    code_to_mic = {i + 1: c for i, c in enumerate(mi.stan_coords.mics)}
    code_to_rxn = {i + 1: c for i, c in enumerate(mi.stan_coords.reactions)}
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
    for var, target_type in zip(
        ["yenz_sim", "yconc_sim", "yflux_sim"], ["enzyme", "mic", "flux"]
    ):
        if var in sim.stan_variables().keys():
            for i, y in enumerate(sim.stan_variable(var)[0]):
                exp_id, var_id = exp_ids[var][i], var_ids[var][i]
                measurement = [
                    m
                    for m in mi.measurements
                    if m.experiment_id == exp_id
                    and m.target_type == target_type
                    and m.target_id == var_id
                ][0]
                measurement.value = y
    return out


def enrich_true_values(tvin, input_data):
    """Add true values for auxiliary parameters."""

    def logz_for_vec(truth, priors):
        if len(truth) == 0:
            return []
        t = np.array(truth)
        lc = np.array(priors[0])
        s = np.array(priors[1])
        return (np.log(t) - np.log(lc)) / s

    def z_for_vec(truth, priors):
        if len(truth) == 0:
            return []
        t = np.array(truth)
        lc = np.array(priors[0])
        s = np.array(priors[1])
        return (t - lc) / s

    def logz_for_mat(truth, priors):
        if len(truth) == 0:
            return np.array([[] for ex in priors])
        out = []
        for i, ex_priors in enumerate(priors):
            t = np.array(truth[i])
            lc = np.array(ex_priors[0])
            s = np.array(ex_priors[1])
            out.append((np.log(t) - np.log(lc)) / s)
        return np.array(out)

    def z_for_mat(truth, priors):
        if len(truth) == 0:
            return np.array([[] for ex in priors])
        out = []
        for i, ex_priors in enumerate(priors):
            t = np.array(truth[i])
            lc = np.array(ex_priors[0])
            s = np.array(ex_priors[1])
            out.append((t - lc) / s)
        return np.array(out)

    return {
        **tvin,
        **{
            "log_km_z": logz_for_vec(tvin["km"], input_data["km_priors"]),
            "log_kcat_z": logz_for_vec(tvin["kcat"], input_data["kcat_priors"]),
            "log_ki_z": logz_for_vec(tvin["ki"], input_data["ki_priors"]),
            "log_dissociation_constant_t_z": logz_for_vec(
                tvin["diss_t"], input_data["diss_t_priors"]
            ),
            "log_dissociation_constant_r_z": logz_for_vec(
                tvin["diss_r"], input_data["diss_r_priors"]
            ),
            "log_transfer_constant_z": logz_for_vec(
                tvin["transfer_constant"], input_data["tc_priors"]
            ),
            "log_enzyme_z": logz_for_mat(tvin["enzyme"], input_data["enzyme_priors"]),
            "log_conc_unbalanced_z": logz_for_mat(
                tvin["conc_unbalanced"], input_data["unbalanced_priors"]
            ),
            "drain_z": logz_for_mat(tvin["drain"], input_data["drain_priors"]),
            "log_phos_kcat_z": logz_for_vec(
                tvin["phos_kcat"], input_data["phos_kcat_priors"]
            ),
            "log_phos_conc_z": logz_for_vec(
                tvin["phos_enzyme_conc"], input_data["phos_conc_priors"]
            ),
        },
    }


def run_simulation_study(mi_in: MaudInput, true_params_raw):
    """Run a simulation study.

    :param mi_in: A MaudInput object
    :param true_params_raw: dictionary of param name -> true param values

    """

    # compile stan model
    here = os.path.dirname(os.path.abspath(__file__))
    stan_path = os.path.join(here, STAN_PROGRAM_RELATIVE_PATH)
    model = CmdStanModel(stan_file=stan_path)
    # generate input data for simulation
    input_data_sim = get_input_data(mi_in)
    # get all true values (including for non-centered parameters)
    true_params = enrich_true_values(true_params_raw, input_data_sim)
    # generate simulated measurements
    sim = model.sample(data=input_data_sim, inits=true_params, **SIM_CONFIG)
    # extract simulated measurements and add them to mi_in
    mi = add_measurements_to_maud_input(mi_in, sim, input_data_sim)
    # create new input data
    input_data_sample = get_input_data(mi)
    # sample
    samples = model.sample(
        data=input_data_sample, inits=true_params, **mi.config.cmdstanpy_config
    )
    return SimulationStudyOutput(
        input_data_sim, input_data_sample, true_params, sim, mi, samples
    )
