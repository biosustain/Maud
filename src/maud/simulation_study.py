"""Code for running simulation studies."""


from copy import deepcopy

import numpy as np
from cmdstanpy import CmdStanMCMC, CmdStanModel

from maud.data_model import MaudInput, SimulationStudyOutput
from maud.sampling import get_input_data


ODE_CONFIG = {
    "abs_tol": 1e-7,
    "rel_tol": 1e-7,
    "max_num_steps": int(1e9),
    "likelihood": 1,
    "timepoint": 500,
}
SIM_CONFIG = {
    "chains": 1,
    "fixed_param": True,
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
                (np.log(np.array(truth[i])) - np.log(np.array(loc[i])))
                / np.array(scale[i])
                for i, _ in enumerate(truth)
            ]
        )

    def z_for_mat(truth, loc, scale):
        return np.array(
            [
                (np.array(truth[i]) - np.array(loc[i])) / np.array(scale[i])
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
                tvin["enzyme"],
                input_data["prior_loc_enzyme"],
                input_data["prior_scale_enzyme"],
            ),
            "log_conc_unbalanced_z": logz_for_mat(
                tvin["conc_unbalanced"],
                input_data["prior_loc_unbalanced"],
                input_data["prior_scale_unbalanced"],
            ),
            "log_drain_z": z_for_mat(
                tvin["drain"],
                input_data["prior_loc_drain"],
                input_data["prior_scale_drain"],
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


def run_simulation_study(
    stan_path: str, mi_in: MaudInput, true_values_in: dict, sample_config: dict
):
    """Run a simulation study."""

    # compile stan model
    model = CmdStanModel(stan_file=stan_path)
    # generate input data for simulation
    input_data_sim = get_input_data(mi_in, **ODE_CONFIG)
    # get all true values (including for non-centered parameters)
    true_values = enrich_true_values(true_values_in, input_data_sim)
    # generate simulated measurements
    sim = model.sample(data=input_data_sim, inits=true_values, **SIM_CONFIG)
    # extract simulated measurements and add them to mi_in
    mi = add_measurements_to_maud_input(mi_in, sim, input_data_sim)
    # create new input data
    input_data_sample = get_input_data(mi, **ODE_CONFIG)
    # sample
    samples = model.sample(data=input_data_sample, inits=true_values, **sample_config)
    return SimulationStudyOutput(
        input_data_sim, input_data_sample, true_values, sim, mi, samples
    )
