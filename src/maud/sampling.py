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
    enzyme_codes = mi.stan_codes.enzyme_codes
    mic_codes = mi.stan_codes.mic_codes
    drain_codes = mi.stan_codes.drain_codes
    reaction_codes = mi.stan_codes.reaction_codes
    S_enz = pd.DataFrame(0, index=enzyme_codes, columns=mic_codes)
    S_drain = pd.DataFrame(0, index=drain_codes, columns=mic_codes)
    S_reactions = pd.DataFrame(0, index=mic_codes, columns=reaction_codes + drain_codes)
    S_enz_to_flux_map = pd.DataFrame(0, index=reaction_codes, columns=enzyme_codes)
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
    balanced_mic_codes = mi.stan_codes.balanced_mic_codes
    rxn_measurements = [m for m in mi.measurements if m.target_type == "flux"]
    complete_reactions = mi.stan_codes.reaction_codes + mi.stan_codes.drain_codes
    for exp in mi.stan_codes.experiment_codes:
        measured_rxns = [
            m.target_id for m in rxn_measurements if m.experiment_id == exp
        ]
        flux_paths = get_null_space(rxn_stoic.loc[balanced_mic_codes].values)
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
    experiment_codes = mi.stan_codes.experiment_codes
    col_codes = (
        mi.stan_codes.enzyme_codes
        if knockout_type == "enz"
        else mi.stan_codes.phos_enz_codes
    )
    knockout_matrix = np.zeros([len(experiment_codes), len(col_codes)])
    for knockout in mi.knockouts:
        if knockout.knockout_type == knockout_type:
            i = codify(experiment_codes)[knockout.experiment_id] - 1
            j = codify(col_codes)[knockout.target_id] - 1
            knockout_matrix[i, j] = 1
    return knockout_matrix.tolist()


def get_phos_act_inh_matrix(mi: MaudInput):
    """Get two enzyme/phos_enzyme matrix, 1 if enyzme modified, else 0.

    :param mi: a MaudInput object
    """
    enzyme_codes = mi.stan_codes.enzyme_codes
    phos_enz_codes = mi.stan_codes.phos_enz_codes
    S_phos_act = np.zeros([len(enzyme_codes), len(phos_enz_codes)])
    S_phos_inh = np.zeros([len(enzyme_codes), len(phos_enz_codes)])
    for phos in mi.kinetic_model.phosphorylation:
        i = codify(enzyme_codes)[phos.enzyme_id] - 1
        j = codify(phos_enz_codes)[phos.id] - 1
        if phos.activating:
            S_phos_act[i, j] = 1
        elif phos.inhbiting:
            S_phos_inh[i, j] = 1
    return S_phos_act.T.tolist(), S_phos_inh.T.tolist()


def _get_km_lookup(mi: MaudInput) -> List[List[int]]:
    out = []
    for mic in mi.stan_codes.mic_codes:
        row = []
        for enz in mi.stan_codes.enzyme_codes:
            km_id = f"{mic}_{enz}"
            if km_id in mi.stan_codes.km_codes:
                row.append(codify(mi.stan_codes.km_codes)[km_id])
            else:
                row.append(0)
        out.append(row)
    return out


def _tabulate_priors_1d(priors: List[Prior]) -> List[List[float]]:
    if len(priors) == 0:
        return [[], []]
    return [[p.location for p in priors], [p.scale for p in priors]]


def _tabulate_priors_2d(priors: List[Prior], exp_codes, target_codes, defaults):
    out = []
    for exp in exp_codes:
        locs = []
        scales = []
        for target in target_codes:
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
    return [
        [0.01 for mic in mi.kinetic_model.mics] for _ in mi.stan_codes.experiment_codes
    ]


def get_input_data(mi: MaudInput) -> dict:
    """Get the input to inference_model.stan from a MaudInput object.

    :param mi: a MaudInput object

    """
    validate_specified_fluxes(mi)
    sorted_enzymes = sorted(
        [e for r in mi.kinetic_model.reactions for e in r.enzymes],
        key=lambda e: codify(mi.stan_codes.enzyme_codes)[e.id],
    )
    water_stoichiometry = [r.water_stoichiometry for r in mi.kinetic_model.reactions]
    mic_to_met = [
        codify(mi.stan_codes.metabolite_codes)[mic.metabolite_id]
        for mic in mi.kinetic_model.mics
    ]
    S_enz, S_to_flux, S_full, S_drain, _ = get_stoics(mi)
    knockout_matrix_enzyme = get_knockout_matrix(mi, knockout_type="enz")
    knockout_matrix_phos = get_knockout_matrix(mi, knockout_type="phos")
    S_phos_act, S_phos_inh = get_phos_act_inh_matrix(mi)
    unbalanced_mic_ix, balanced_mic_ix = (
        [codify(mi.stan_codes.mic_codes)[mic_id] for mic_id in codes]
        for codes in (
            mi.stan_codes.unbalanced_mic_codes,
            mi.stan_codes.balanced_mic_codes,
        )
    )

    return {
        # sizes
        "N_mic": len(mi.kinetic_model.mics),
        "N_unbalanced": len(mi.stan_codes.unbalanced_mic_codes),
        "N_metabolite": len(mi.stan_codes.metabolite_codes),
        "N_km": len(mi.stan_codes.km_codes),
        "N_reaction": len(mi.stan_codes.reaction_codes),
        "N_enzyme": len(mi.stan_codes.enzyme_codes),
        "N_phosphorylation_enzymes": len(mi.stan_codes.phos_enz_codes),
        "N_experiment": len(mi.stan_codes.experiment_codes),
        "N_flux_measurement": len(mi.stan_codes.yflux_reaction_codes),
        "N_enzyme_measurement": len(mi.stan_codes.yenz_enz_codes),
        "N_conc_measurement": len(mi.stan_codes.yconc_mic_codes),
        "N_ki": len(mi.priors.inhibition_constant_priors),
        "N_ai": len(mi.priors.tense_dissociation_constant_priors),
        "N_aa": len(mi.priors.relaxed_dissociation_constant_priors),
        "N_ae": len(mi.priors.transfer_constant_priors),
        "N_drain": len(mi.stan_codes.drain_codes),
        # codes
        "unbalanced_mic_ix": unbalanced_mic_ix,
        "balanced_mic_ix": balanced_mic_ix,
        "experiment_yconc": mi.stan_codes.yconc_exp_codes,
        "mic_ix_yconc": mi.stan_codes.yconc_mic_codes,
        "experiment_yflux": mi.stan_codes.yflux_exp_codes,
        "reaction_yflux": mi.stan_codes.yflux_reaction_codes,
        "experiment_yenz": mi.stan_codes.yenz_exp_codes,
        "enzyme_yenz": mi.stan_codes.yenz_enz_codes,
        "ci_ix": mi.stan_codes.ci_mic_codes,
        "ai_ix": mi.stan_codes.ai_mic_codes,
        "aa_ix": mi.stan_codes.aa_mic_codes,
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
            mi.stan_codes.experiment_codes,
            mi.stan_codes.unbalanced_mic_codes,
            [DEFAULT_PRIOR_LOC_UNBALANCED, DEFAULT_PRIOR_SCALE_UNBALANCED],
        ),
        "enzyme_priors": _tabulate_priors_2d(
            mi.priors.enzyme_concentration_priors,
            mi.stan_codes.experiment_codes,
            mi.stan_codes.enzyme_codes,
            [DEFAULT_PRIOR_LOC_ENZYME, DEFAULT_PRIOR_SCALE_ENZYME],
        ),
        "phos_conc_priors": _tabulate_priors_2d(
            mi.priors.phos_enz_concentration_priors,
            mi.stan_codes.experiment_codes,
            mi.stan_codes.phos_enz_codes,
            [DEFAULT_PRIOR_LOC_ENZYME, DEFAULT_PRIOR_SCALE_ENZYME],
        ),
        "drain_priors": _tabulate_priors_2d(
            mi.priors.drain_priors,
            mi.stan_codes.experiment_codes,
            mi.stan_codes.drain_codes,
            [0, 1],
        ),
        # config
        "rel_tol": mi.config.ode_config["rel_tol"],
        "abs_tol": mi.config.ode_config["abs_tol"],
        "max_num_steps": int(mi.config.ode_config["max_num_steps"]),
        "LIKELIHOOD": int(mi.config.likelihood),
        "timepoint": mi.config.ode_config["timepoint"],
        "conc_init": _get_conc_init(mi),
    }
