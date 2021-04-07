# Copyright (C) 2019 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Code for sampling from a posterior distribution."""

import os
import warnings
from typing import List

import cmdstanpy
import numpy as np
import pandas as pd

from maud.data_model import MaudInput, Prior
from maud.utils import codify, get_null_space, get_rref


HERE = os.path.dirname(os.path.abspath(__file__))
INCLUDE_PATH = ""
DEFAULT_PRIOR_LOC_UNBALANCED = 0.1
DEFAULT_PRIOR_SCALE_UNBALANCED = 2.0
DEFAULT_PRIOR_LOC_ENZYME = 0.1
DEFAULT_PRIOR_SCALE_ENZYME = 2.0
DEFAULT_PRIOR_LOC_DRAIN = None
DEFAULT_PRIOR_SCALE_DRAIN = None
STAN_PROGRAM_RELATIVE_PATH = "inference_model.stan"

DEFAULT_SAMPLE_CONFIG = {
    "iter_warmup": 5,
    "iter_sampling": 5,
    "chains": 2,
    "max_treedepth": 11,
    "inits": 0,
    "show_progress": True,
    "step_size": 0.025,
    "adapt_delta": 0.99,
    "save_warmup": True,
    "threads_per_chain": 1,
}
DEFAULT_ODE_CONFIG = {
    "rel_tol": 1e-9,
    "abs_tol": 1e-9,
    "max_num_steps": int(1e9),
    "timepoint": 500,
}
SIM_CONFIG = {
    "chains": 1,
    "fixed_param": True,
    "inits": 0,
    "iter_warmup": 0,
    "show_progress": False,
    "threads_per_chain": 1,
}


def sample(mi: MaudInput, output_dir: str) -> cmdstanpy.CmdStanMCMC:
    """Sample from the posterior defined by mi.

    :param mi: a MaudInput object
    :param output_dir: a string specifying where to save the output.
    """
    config = {
        **DEFAULT_SAMPLE_CONFIG,
        **mi.config.cmdstanpy_config,
        **{"output_dir": output_dir},
    }
    return _sample_given_config(mi, output_dir, config)


def simulate(mi: MaudInput, output_dir: str, n: int) -> cmdstanpy.CmdStanMCMC:
    """Generate simulations from the prior mean.

    :param mi: a MaudInput object
    :param output_dir: a string specifying where to save the output.
    """
    config = {**SIM_CONFIG, **{"output_dir": output_dir, "iter_sampling": n}}
    return _sample_given_config(mi, output_dir, config)


def _sample_given_config(
    mi: MaudInput, output_dir: str, config: dict
) -> cmdstanpy.CmdStanMCMC:
    """Call CmdStanModel.sample, having already specified all arguments.

    :param mi: a MaudInput object
    :param output_dir: a string specifying where to save the output.
    :param config: a dictionary of keyword arguments to CmdStanModel.sample.
    """

    input_filepath = os.path.join(output_dir, "input_data.json")
    input_data = get_input_data(mi)
    cmdstanpy.utils.jsondump(input_filepath, input_data)
    stan_program_filepath = os.path.join(HERE, STAN_PROGRAM_RELATIVE_PATH)
    include_path = os.path.join(HERE, INCLUDE_PATH)
    cpp_options = {}
    stanc_options = {"include_paths": [include_path]}
    if config["threads_per_chain"] != 1:
        cpp_options["STAN_THREADS"] = True
        os.environ["STAN_NUM_THREADS"] = str(config["threads_per_chain"])
    model = cmdstanpy.CmdStanModel(
        stan_file=stan_program_filepath,
        stanc_options=stanc_options,
        cpp_options=cpp_options,
    )
    return model.sample(data=input_filepath, **config)


def get_stoics(mi: MaudInput):
    """Get full stoichiometric matrix for each isoenzyme.

    :param mi: A MaudInput object
    """
    enzyme_ix = mi.stan_coords.enzymes
    mic_ix = mi.stan_coords.mics
    drain_ix = mi.stan_coords.drains
    reaction_ix = mi.stan_coords.reactions
    S_enz = pd.DataFrame(0, index=enzyme_ix, columns=mic_ix)
    S_drain = pd.DataFrame(0, index=drain_ix, columns=mic_ix)
    S_reactions = pd.DataFrame(0, index=mic_ix, columns=reaction_ix + drain_ix)
    S_enz_to_flux_map = pd.DataFrame(0, index=reaction_ix, columns=enzyme_ix)
    for rxn in mi.kinetic_model.reactions:
        for enz in rxn.enzymes:
            for met, stoic in rxn.stoichiometry.items():
                S_enz_to_flux_map.loc[rxn.id, enz.id] = 1
                S_enz.loc[enz.id, met] = stoic
                S_reactions.loc[met, rxn.id] = stoic
    if S_drain.shape[0] > 0:
        for drain in mi.kinetic_model.drains:
            for met, stoic in drain.stoichiometry.items():
                S_drain.loc[drain.id, met] = stoic
                S_reactions.loc[met, drain.id] = stoic
    S_full = pd.concat([S_enz, S_drain], ignore_index=True)
    return S_enz, S_enz_to_flux_map, S_full, S_drain, S_reactions


def validate_specified_fluxes(mi: MaudInput):
    """Check that appropriate fluxes have been measured.

    :param mi: A MaudInput object
    """
    _, _, _, _, rxn_stoic = get_stoics(mi)
    balanced_mic_ix = mi.stan_coords.balanced_mics
    rxn_measurements = [m for m in mi.measurements if m.target_type == "flux"]
    complete_reactions = mi.stan_coords.reactions + mi.stan_coords.drains
    for exp in mi.stan_coords.experiments:
        measured_rxns = [
            m.target_id for m in rxn_measurements if m.experiment_id == exp
        ]
        flux_paths = get_null_space(rxn_stoic.loc[balanced_mic_ix].values)
        _, n_dof = np.shape(flux_paths)
        rref_flux_paths = np.matrix(get_rref(flux_paths.T))
        rref_flux_paths[np.abs(rref_flux_paths) < 1e-10] = 0
        flux_paths_df = pd.DataFrame(rref_flux_paths, columns=complete_reactions)
        for _, flux_path in flux_paths_df.iterrows():
            if any(flux_path[measured_rxns]) != 0:
                pass
            else:
                possible_measurements = []
                for rxn, st in flux_path.items():
                    if st != 0:
                        possible_measurements.append(rxn)
                msg = (
                    "\nYour system appears to be underdetermined in "
                    + f"experiment: {exp}\n"
                    + "Please define a reaction from the following list:\n"
                    + "\n".join(possible_measurements)
                )
                warnings.warn(msg)

        if len(measured_rxns) > n_dof:
            msg = (
                "You appear to have specified too many reactions.\n"
                + "This will bias the statistical model\n"
                + "as the measurements are not independent."
            )
            warnings.warn(msg)


def get_knockout_matrix(mi: MaudInput, knockout_type: str):
    """Get binary experiment, enzyme matrix, 1 if enyzme knocked out, 0 if not.

    :param mi: a MaudInput object
    :param knockout_type: either "enz" or "phos"
    """
    if knockout_type not in ["enz", "phos"]:
        raise ValueError("knockout_type must be either 'enz' or 'phos'.")
    experiment_ix = mi.stan_coords.experiments
    col_ix = (
        mi.stan_coords.enzymes if knockout_type == "enz" else mi.stan_coords.phos_enzs
    )
    knockout_matrix = np.zeros([len(experiment_ix), len(col_ix)])
    for knockout in mi.knockouts:
        if knockout.knockout_type == knockout_type:
            i = codify(experiment_ix)[knockout.experiment_id] - 1
            j = codify(col_ix)[knockout.target_id] - 1
            knockout_matrix[i, j] = 1
    return knockout_matrix.tolist()


def get_phos_act_inh_matrix(mi: MaudInput):
    """Get two enzyme/phos_enzyme matrix, 1 if enyzme modified, else 0.

    :param mi: a MaudInput object
    """
    enzyme_ix = mi.stan_coords.enzymes
    phos_enz_ix = mi.stan_coords.phos_enzs
    S_phos_act = np.zeros([len(enzyme_ix), len(phos_enz_ix)])
    S_phos_inh = np.zeros([len(enzyme_ix), len(phos_enz_ix)])
    for phos in mi.kinetic_model.phosphorylation:
        i = codify(enzyme_ix)[phos.enzyme_id] - 1
        j = codify(phos_enz_ix)[phos.id] - 1
        if phos.activating:
            S_phos_act[i, j] = 1
        elif phos.inhbiting:
            S_phos_inh[i, j] = 1
    return S_phos_act.T.tolist(), S_phos_inh.T.tolist()


def _get_km_lookup(mi: MaudInput) -> List[List[int]]:
    out = []
    for mic in mi.stan_coords.mics:
        row = []
        for enz in mi.stan_coords.enzymes:
            km_id = f"{enz}_{mic}"
            if km_id in mi.stan_coords.kms:
                row.append(codify(mi.stan_coords.kms)[km_id])
            else:
                row.append(0)
        out.append(row)
    return out


def _tabulate_priors_1d(priors: List[Prior]) -> List[List[float]]:
    if len(priors) == 0:
        return [[], []]
    return [[p.location for p in priors], [p.scale for p in priors]]


def _tabulate_priors_2d(priors: List[Prior], exp_ix, target_ix, defaults):
    out = []
    for exp in exp_ix:
        locs = []
        scales = []
        for target in target_ix:
            match = [
                p
                for p in priors
                if p.experiment_id == exp and target in p.__dict__.values()
            ]
            if any(match):
                p = match[0]
                locs.append(p.location)
                scales.append(p.scale)
            else:
                locs.append(defaults[0])
                scales.append(defaults[1])
        out.append([locs, scales])
    return out


def _get_conc_init(mi):
    """Get the initial concentrations for ODE solver from MaudInput object.

    :param mi: a MaudInput object

    """
    conc_init = [
        [0.01 for mic in mi.kinetic_model.mics] for _ in mi.stan_coords.experiments
    ]
    for row in mi.measurements:
        if row.target_type == "mic":
            if row.target_id in mi.stan_coords.balanced_mics:
                mic_idx = codify(mi.stan_coords.mics)[row.target_id] - 1
                exp_idx = codify(mi.stan_coords.experiments)[row.experiment_id] - 1
                conc_init[exp_idx][mic_idx] = row.value
    for p in mi.priors.unbalanced_metabolite_priors:
        mic_idx = codify(mi.stan_coords.mics)[p.mic_id] - 1
        exp_idx = codify(mi.stan_coords.experiments)[p.experiment_id] - 1
        conc_init[exp_idx][mic_idx] = p.location
    return conc_init


def get_input_data(mi: MaudInput) -> dict:
    """Get the input to inference_model.stan from a MaudInput object.

    :param mi: a MaudInput object

    """
    validate_specified_fluxes(mi)
    sorted_enzymes = sorted(
        [e for r in mi.kinetic_model.reactions for e in r.enzymes],
        key=lambda e: codify(mi.stan_coords.enzymes)[e.id],
    )
    water_stoichiometry = [e.water_stoichiometry for e in sorted_enzymes]
    mic_to_met = [
        codify(mi.stan_coords.metabolites)[mic.metabolite_id]
        for mic in mi.kinetic_model.mics
    ]
    S_enz, S_to_flux, S_full, S_drain, _ = get_stoics(mi)
    knockout_matrix_enzyme = get_knockout_matrix(mi, knockout_type="enz")
    knockout_matrix_phos = get_knockout_matrix(mi, knockout_type="phos")
    S_phos_act, S_phos_inh = get_phos_act_inh_matrix(mi)
    unbalanced_mic_ix, balanced_mic_ix = (
        [codify(mi.stan_coords.mics)[mic_id] for mic_id in ix]
        for ix in (
            mi.stan_coords.unbalanced_mics,
            mi.stan_coords.balanced_mics,
        )
    )
    allosteric_enzymes = [
        e
        for e in sorted_enzymes
        if (
            len(e.modifiers["allosteric_inhibitor"]) > 0
            or len(e.modifiers["allosteric_activator"]) > 0
        )
    ]
    exp_to_ix = codify(mi.stan_coords.experiments)
    mic_to_ix = codify(mi.stan_coords.mics)
    rxn_to_ix = codify(mi.stan_coords.reactions)
    enz_to_ix = codify(mi.stan_coords.enzymes)
    return {
        # sizes
        "N_mic": len(mi.kinetic_model.mics),
        "N_unbalanced": len(mi.stan_coords.unbalanced_mics),
        "N_metabolite": len(mi.stan_coords.metabolites),
        "N_km": len(mi.stan_coords.kms),
        "N_reaction": len(mi.stan_coords.reactions),
        "N_enzyme": len(mi.stan_coords.enzymes),
        "N_phosphorylation_enzymes": len(mi.stan_coords.phos_enzs),
        "N_experiment": len(mi.stan_coords.experiments),
        "N_flux_measurement": len(mi.stan_coords.yflux_rxns),
        "N_enzyme_measurement": len(mi.stan_coords.yenz_enzs),
        "N_conc_measurement": len(mi.stan_coords.yconc_mics),
        "N_ki": len(mi.stan_coords.ci_mics),
        "N_ai": len(mi.stan_coords.ai_mics),
        "N_aa": len(mi.stan_coords.aa_mics),
        "N_ae": len(allosteric_enzymes),
        "N_drain": len(mi.stan_coords.drains),
        # indexes
        "unbalanced_mic_ix": unbalanced_mic_ix,
        "balanced_mic_ix": balanced_mic_ix,
        "experiment_yconc": [exp_to_ix[e] for e in mi.stan_coords.yconc_exps],
        "mic_ix_yconc": [mic_to_ix[m] for m in mi.stan_coords.yconc_mics],
        "experiment_yflux": [exp_to_ix[e] for e in mi.stan_coords.yflux_exps],
        "reaction_yflux": [rxn_to_ix[r] for r in mi.stan_coords.yflux_rxns],
        "experiment_yenz": [exp_to_ix[e] for e in mi.stan_coords.yenz_exps],
        "enzyme_yenz": [enz_to_ix[e] for e in mi.stan_coords.yenz_enzs],
        "ci_ix": [mic_to_ix[m] for m in mi.stan_coords.ci_mics],
        "ai_ix": [mic_to_ix[m] for m in mi.stan_coords.ai_mics],
        "aa_ix": [mic_to_ix[m] for m in mi.stan_coords.aa_mics],
        # network properties
        "S_enz": S_enz.T.values,
        "S_to_flux_map": S_to_flux.values,
        "S_drain": S_drain.T.values,
        "S_full": S_full.T.values,
        "S_phos_act": S_phos_act,
        "S_phos_inh": S_phos_inh,
        "water_stoichiometry": water_stoichiometry,
        "mic_to_met": mic_to_met,
        "km_lookup": _get_km_lookup(mi),
        "is_knockout": knockout_matrix_enzyme,
        "is_phos_knockout": knockout_matrix_phos,
        "subunits": [e.subunits for e in sorted_enzymes],
        "n_ci": [len(e.modifiers["competitive_inhibitor"]) for e in sorted_enzymes],
        "n_ai": [len(e.modifiers["allosteric_inhibitor"]) for e in sorted_enzymes],
        "n_aa": [len(e.modifiers["allosteric_activator"]) for e in sorted_enzymes],
        # measurements
        "yconc": [m.value for m in mi.measurements if m.target_type == "mic"],
        "sigma_conc": [m.error for m in mi.measurements if m.target_type == "mic"],
        "yflux": [m.value for m in mi.measurements if m.target_type == "flux"],
        "sigma_flux": [m.error for m in mi.measurements if m.target_type == "flux"],
        "yenz": [m.value for m in mi.measurements if m.target_type == "enzyme"],
        "sigma_enz": [m.error for m in mi.measurements if m.target_type == "enzyme"],
        # priors
        "fe_priors": _tabulate_priors_1d(mi.priors.formation_energy_priors),
        "kcat_priors": _tabulate_priors_1d(mi.priors.kcat_priors),
        "km_priors": _tabulate_priors_1d(mi.priors.km_priors),
        "ki_priors": _tabulate_priors_1d(mi.priors.inhibition_constant_priors),
        "diss_t_priors": _tabulate_priors_1d(
            mi.priors.tense_dissociation_constant_priors
        ),
        "diss_r_priors": _tabulate_priors_1d(
            mi.priors.relaxed_dissociation_constant_priors
        ),
        "phos_kcat_priors": _tabulate_priors_1d(mi.priors.phos_kcat_priors),
        "tc_priors": _tabulate_priors_1d(mi.priors.transfer_constant_priors),
        "unbalanced_priors": _tabulate_priors_2d(
            mi.priors.unbalanced_metabolite_priors,
            mi.stan_coords.experiments,
            mi.stan_coords.unbalanced_mics,
            [DEFAULT_PRIOR_LOC_UNBALANCED, DEFAULT_PRIOR_SCALE_UNBALANCED],
        ),
        "enzyme_priors": _tabulate_priors_2d(
            mi.priors.enzyme_concentration_priors,
            mi.stan_coords.experiments,
            mi.stan_coords.enzymes,
            [DEFAULT_PRIOR_LOC_ENZYME, DEFAULT_PRIOR_SCALE_ENZYME],
        ),
        "phos_conc_priors": _tabulate_priors_2d(
            mi.priors.phos_enz_concentration_priors,
            mi.stan_coords.experiments,
            mi.stan_coords.phos_enzs,
            [DEFAULT_PRIOR_LOC_ENZYME, DEFAULT_PRIOR_SCALE_ENZYME],
        ),
        "drain_priors": _tabulate_priors_2d(
            mi.priors.drain_priors,
            mi.stan_coords.experiments,
            mi.stan_coords.drains,
            [DEFAULT_PRIOR_LOC_DRAIN, DEFAULT_PRIOR_SCALE_DRAIN],
        ),
        # config
        "rel_tol": mi.config.ode_config["rel_tol"],
        "abs_tol": mi.config.ode_config["abs_tol"],
        "max_num_steps": int(mi.config.ode_config["max_num_steps"]),
        "LIKELIHOOD": int(mi.config.likelihood),
        "timepoint": mi.config.ode_config["timepoint"],
        "conc_init": _get_conc_init(mi),
    }
