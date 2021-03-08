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
from typing import Dict, List

import cmdstanpy
import numpy as np
import pandas as pd

from maud.data_model import KineticModel, MaudInput, Prior
from maud.utils import get_null_space, get_rref, codify



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

    :param kinetic_model: A Kinetic Model object
    :param enzyme_codes: the codified enzyme codes
    :param metabolite_codes: the codified metabolite codes
    """
    kinetic_model = mi.kinetic_model
    enzyme_codes = mi.stan_codes.enzyme_codes
    metabolite_codes = mi.stan_codes.metabolite_codes
    mic_codes = mi.stan_codes.mic_codes
    reaction_codes = mi.stan_codes.reaction_codes
    drain_codes = mi.stan_codes.drain_codes
    S_enz = pd.DataFrame(0, index=enzyme_codes, columns=mic_codes)
    S_reactions = pd.DataFrame(
        0, index=mic_codes, columns=reaction_codes+drain_codes
    )
    if len(drain_codes) > 0:
        S_drain = pd.DataFrame(0, index=drain_codes, columns=mic_codes)
        S_complete = pd.DataFrame(
            0, index=enzyme_codes+drain_codes, columns=mic_codes
        )
    else:
        S_complete = pd.DataFrame(0, index=enzyme_codes, columns=mic_codes)
    S_enz_to_flux_map = pd.DataFrame(0, index=reaction_codes, columns=enzyme_codes)
    for rxn in kinetic_model.reactions:
        for enz in rxn.enzymes:
            for met, stoic in rxn.stoichiometry.items():
                S_enz_to_flux_map.loc[rxn.id, enz.id] = 1
                S_enz.loc[enz.id, met] = stoic
                S_complete.loc[enz.id, met] = stoic
                S_reactions.loc[met, rxn.id] = stoic
    if len(drain_codes) > 0:
        for drain in kinetic_model.drains:
            for met, stoic in drain.stoichiometry.items():
                S_drain.loc[drain.id, met] = stoic
                S_complete.loc[drain.id, met] = stoic
                S_reactions.loc[met, drain.id] = stoic
    else:
        S_drain = pd.DataFrame()
    return S_enz, S_enz_to_flux_map, S_complete, S_drain, S_reactions


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
        n_rxns, n_dof = np.shape(flux_paths)
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
        j = codify(phos_enz_codes)[phos_id] - 1
        if phos.activating:
            S_phos_act[i, j] = 1
        elif phos.inhbiting:
            S_phos_inh[i, j] = 1
    return S_phos_act.tolist(), S_phos_inh.tolist()


def get_km_lookup(mi: MaudInput) -> List[List[int]]:
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


def tabulate_priors_1d(priors: List[Prior]) -> List[List[float]]:
    if len(priors) == 0:
        return [[]]
    return [[p.location, p.scale] for p in priors]

def tabulate_priors_2d(priors: List[Prior], exp_codes, target_codes, defaults):
    out = []
    for exp in exp_codes:
        row = []
        for target in target_codes:
            f = lambda p: p.experiment_id == exp and p.target_id == target
            match = list(filter(f, priors))
            if any(match):
                p = match[0]
                row.append([p.location, p.scale])
            else:
                row.append(defaults)
        out.append(row)
    return out


def get_conc_init(mi):
    return [
        [0.01 for mic in mi.kinetic_model.mics]
        for _ in mi.stan_codes.experiment_codes
    ]


def get_input_data(mi: MaudInput) -> dict:
    """Re-write of this function under the following assumptions:
    
    - default priors are taken care of at the io stage, so mi.priors is
      complete.
    - mi.stan_codes includes yenz_enz_codes, yconc_mic_codes, etc.
    - stan codes are lists
    - mi has an attribut measurements 

    """
    validate_specified_fluxes(mi)
    sorted_enzymes = sorted(
        [e for r in mi.kinetic_model.reactions for e in r.enzymes],
        key=lambda e: codify(mi.stan_codes.enzyme_codes)[e.id]
    )
    water_stoichiometry = [
        r.water_stoichiometry for r in mi.kinetic_model.reactions
    ]
    mic_to_met = [
        codify(mi.stan_codes.metabolite_codes)[mic.metabolite_id]
        for mic in mi.kinetic_model.mics
    ]
    S_enz, S_to_flux, S_full, S_drain, S_reactions = get_stoics(mi)
    knockout_matrix_enzyme = get_knockout_matrix(mi, knockout_type="enz")
    knockout_matrix_phos = get_knockout_matrix(mi, knockout_type="phos")
    S_phos_act, S_phos_inh = get_phos_act_inh_matrix(mi)
    unbalanced_mic_ix, balanced_mic_ix = (
        [codify(mi.stan_codes.mic_codes)[mic_id] for mic_id in codes]
        for codes in (
            mi.stan_codes.unbalanced_mic_codes, mi.stan_codes.balanced_mic_codes
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
        "S_drain": S_drain.values,
        "S_full": S_full.T.values,
        "S_phos_act": S_phos_act,
        "S_phos_inh": S_phos_inh,
        "water_stoichiometry": water_stoichiometry,
        "mic_to_met": mic_to_met,
        "km_lookup": get_km_lookup(mi),
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
        "fe_priors": tabulate_priors_1d(mi.priors.formation_energy_priors),
        "kcat_priors": tabulate_priors_1d(mi.priors.kcat_priors),
        "km_priors": tabulate_priors_1d(mi.priors.km_priors),
        "ki_priors": tabulate_priors_1d(mi.priors.inhibition_constant_priors),
        "diss_t_priors": tabulate_priors_1d(mi.priors.tense_dissociation_constant_priors),
        "diss_r_priors": tabulate_priors_1d(mi.priors.relaxed_dissociation_constant_priors),
        "phos_kcat_priors": tabulate_priors_1d(mi.priors.phos_kcat_priors),
        "tc_priors": tabulate_priors_1d(mi.priors.transfer_constant_priors),
        "unbalanced_priors": tabulate_priors_2d(
            mi.priors.unbalanced_metabolite_priors,
            mi.stan_codes.experiment_codes,
            mi.stan_codes.unbalanced_mic_codes,
            [DEFAULT_PRIOR_LOC_UNBALANCED, DEFAULT_PRIOR_SCALE_UNBALANCED],
        ),
        "enzyme_priors": tabulate_priors_2d(
            mi.priors.enzyme_concentration_priors,
            mi.stan_codes.experiment_codes,
            mi.stan_codes.enzyme_codes,
            [DEFAULT_PRIOR_LOC_ENZYME, DEFAULT_PRIOR_SCALE_ENZYME],
        ),
        "phos_conc_priors": tabulate_priors_2d(
            mi.priors.phos_enz_concentration_priors,
            mi.stan_codes.experiment_codes,
            mi.stan_codes.phos_enz_codes,
            [DEFAULT_PRIOR_LOC_ENZYME, DEFAULT_PRIOR_SCALE_ENZYME],
        ),
        "drain_priors": tabulate_priors_2d(
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
        "conc_init": get_conc_init(mi),
    }
    



# def get_input_data(mi: MaudInput) -> dict:
#     """Put a MaudInput into a Stan-friendly dictionary.

#     :param mi: a MaudInput object
#     """
#     likelihood = mi.config.likelihood
#     ode_config = {**DEFAULT_ODE_CONFIG, **mi.config.ode_config}

#     metabolites = mi.kinetic_model.metabolites
#     mics = mi.kinetic_model.mics
#     reactions = mi.kinetic_model.reactions
#     reaction_codes = mi.stan_codes.reaction_codes
#     drain_codes = mi.stan_codes.drain_codes
#     enzyme_codes = mi.stan_codes.enzyme_codes
#     experiment_codes = mi.stan_codes.experiment_codes
#     met_codes = mi.stan_codes.metabolite_codes
#     mic_codes = mi.stan_codes.mic_codes
#     phos_enz_codes = mi.stan_codes.phos_enz_codes
#     balanced_mic_codes = mi.stan_codes.balanced_mic_codes
#     unbalanced_mic_codes = mi.stan_codes.unbalanced_mic_codes
#     mic_to_met = {
#         mic_codes[mic.id]: met_codes[mic.metabolite_id] for mic in mics.values()
#     }
#     enzymes = {k: v for r in reactions.values() for k, v in r.enzymes.items()}
#     water_stoichiometry = [
#         r.water_stoichiometry for r in reactions.values() for e in r.enzymes.items()
#     ]
#     (
#         enzyme_stoic,
#         stoic_map_to_flux,
#         full_stoic,
#         drain_stoic,
#         rxn_stoic,
#     ) = get_full_stoichiometry(
#         mi.kinetic_model, enzyme_codes, mic_codes, reaction_codes, drain_codes
#     )
#     subunits = (
#         pd.DataFrame(
#             {
#                 "enzyme_id": [e.id for e in enzymes.values()],
#                 "subunits": [e.subunits for e in enzymes.values()],
#                 "index": [enzyme_codes[e.id] for e in enzymes.values()],
#             }
#         )
#         .sort_values(by="index")
#         .reset_index(drop=True)
#     )

#     # priors
#     prior_dfs = {
#         prior_type: pd.DataFrame(map(lambda p: p.__dict__, priors))
#         if len(priors) > 0
#         else pd.DataFrame([], columns=["location", "scale", "mic_code"])
#         for prior_type, priors in mi.priors.__dict__.items()
#     }
#     prior_dfs["formation_energy_priors"] = (
#         prior_dfs["formation_energy_priors"]
#         .assign(stan_code=lambda df: df["metabolite_id"].map(met_codes))
#         .sort_values("stan_code")
#     )
#     if len(prior_dfs["phos_kcat_priors"]) > 0:
#         prior_dfs["phos_kcat_priors"] = (
#             prior_dfs["phos_kcat_priors"]
#             .assign(stan_code=lambda df: df["phos_enz_id"].map(phos_enz_codes))
#             .sort_values("stan_code")
#         )
#     n_modifier = {}
#     for df_name, param_type in zip(
#         [
#             "inhibition_constant_priors",
#             "tense_dissociation_constant_priors",
#             "relaxed_dissociation_constant_priors",
#         ],
#         ["ki", "diss_t", "diss_r"],
#     ):
#         if len(prior_dfs[df_name]) > 0:
#             df = prior_dfs[df_name]
#             df["mic_code"] = df["mic_id"].map(mic_codes)
#             df["enzyme_code"] = df["enzyme_id"].map(enzyme_codes)
#             df.sort_values("enzyme_code", inplace=True)
#             n_modifier[param_type] = (
#                 df.groupby("enzyme_code")
#                 .size()
#                 .reindex(enzyme_codes.values())
#                 .fillna(0)
#                 .astype(int)
#                 .tolist()
#             )
#         else:
#             n_modifier[param_type] = [0] * len(enzymes)
#     if mi.kinetic_model.drains is not None:
#         prior_loc_drain, prior_scale_drain = (
#             prior_dfs["drain_priors"]
#             .set_index(["experiment_id", "drain_id"])[col]
#             .unstack()
#             .loc[experiment_codes.keys(), drain_codes.keys()]
#             for col in ["location", "scale"]
#         )
#     else:
#         prior_loc_drain = prior_scale_drain = pd.DataFrame([])
#     prior_loc_unb = pd.DataFrame(
#         DEFAULT_PRIOR_LOC_UNBALANCED,
#         index=experiment_codes,
#         columns=unbalanced_mic_codes,
#     )
#     prior_scale_unb = pd.DataFrame(
#         DEFAULT_PRIOR_SCALE_UNBALANCED,
#         index=experiment_codes,
#         columns=unbalanced_mic_codes,
#     )
#     prior_loc_enzyme = pd.DataFrame(
#         DEFAULT_PRIOR_LOC_ENZYME, index=experiment_codes, columns=enzyme_codes
#     )
#     prior_scale_enzyme = pd.DataFrame(
#         DEFAULT_PRIOR_SCALE_ENZYME, index=experiment_codes, columns=enzyme_codes
#     )
#     prior_loc_phos_conc = pd.DataFrame(
#         DEFAULT_PRIOR_LOC_ENZYME, index=experiment_codes, columns=phos_enz_codes
#     )
#     prior_scale_phos_conc = pd.DataFrame(
#         DEFAULT_PRIOR_SCALE_ENZYME, index=experiment_codes, columns=phos_enz_codes
#     )
#     conc_init = pd.DataFrame(
#         0.01, index=experiment_codes.values(), columns=mic_codes.values()
#     )
#     enz_conc_init = pd.DataFrame(
#         DEFAULT_PRIOR_LOC_ENZYME,
#         index=experiment_codes.values(),
#         columns=enzyme_codes.values(),
#     )
#     # measurements
#     mic_measurements, reaction_measurements, enzyme_measurements = (
#         pd.DataFrame(
#             [
#                 [exp.id, meas.target_id, meas.value, meas.uncertainty]
#                 for exp in mi.experiments.experiments
#                 for meas in exp.measurements[measurement_type].values()
#             ],
#             columns=["experiment_id", "target_id", "value", "uncertainty"],
#         )
#         for measurement_type in ["mic", "flux", "enzyme"]
#     )
#     for prior_unb in mi.priors.unbalanced_metabolite_priors:
#         mic_id = mic_codes[prior_unb.mic_id]
#         exp_id = experiment_codes[prior_unb.experiment_id]
#         conc_init.loc[exp_id, mic_id] = prior_unb.location
#     for prior_enz in mi.priors.enzyme_concentration_priors:
#         enz_id = enzyme_codes[prior_enz.enzyme_id]
#         exp_id = experiment_codes[prior_enz.experiment_id]
#         enz_conc_init.loc[exp_id, enz_id] = prior_enz.location
#     for _i, row in mic_measurements.iterrows():
#         if row["target_id"] in balanced_mic_codes.keys():
#             row_ix = experiment_codes[row["experiment_id"]]
#             column_ix = mic_codes[row["target_id"]]
#             conc_init.loc[row_ix, column_ix] = row["value"]
#     for _i, row in enzyme_measurements.iterrows():
#         row_ix = experiment_codes[row["experiment_id"]]
#         column_ix = enzyme_codes[row["target_id"]]
#         enz_conc_init.loc[row_ix, column_ix] = row["value"]
#     knockout_matrix = get_knockout_matrix(mi=mi)
#     phos_knockout_matrix = get_phos_knockout_matrix(mi=mi)
#     if len(phos_enz_codes) > 0:
#         S_phos_act, S_phos_inh = get_phos_act_inh_matrix(mi=mi)
#     else:
#         S_phos_act = pd.DataFrame()
#         S_phos_inh = pd.DataFrame()
#     km_lookup = get_km_lookup(prior_dfs["km_priors"], mic_codes, enzyme_codes)
#     for prior_unb in mi.priors.unbalanced_metabolite_priors:
#         mic_id = prior_unb.mic_id
#         exp_id = prior_unb.experiment_id
#         prior_loc_unb.loc[exp_id, mic_id] = prior_unb.location
#         prior_scale_unb.loc[exp_id, mic_id] = prior_unb.scale
#     for prior_enz in mi.priors.enzyme_concentration_priors:
#         enz_id = prior_enz.enzyme_id
#         exp_id = prior_enz.experiment_id
#         prior_loc_enzyme.loc[exp_id, enz_id] = prior_enz.location
#         prior_scale_enzyme.loc[exp_id, enz_id] = prior_enz.scale
#     for prior_phos_enz in mi.priors.phos_enz_concentration_priors:
#         penz_id = prior_phos_enz.phos_enz_id
#         exp_id = prior_phos_enz.experiment_id
#         prior_loc_phos_conc.loc[exp_id, penz_id] = prior_phos_enz.location
#         prior_scale_phos_conc.loc[exp_id, penz_id] = prior_phos_enz.scale
#     validate_specified_fluxes(
#         rxn_stoic,
#         reaction_measurements,
#         experiment_codes,
#         reaction_codes,
#         drain_codes,
#         balanced_mic_codes,
#     )
#     return {
#         "N_mic": len(mics),
#         "N_unbalanced": len(unbalanced_mic_codes),
#         "N_metabolite": len(metabolites),
#         "N_km": len(prior_dfs["km_priors"]),
#         "N_reaction": len(reactions),
#         "N_enzyme": len(enzymes),
#         "N_phosphorylation_enzymes": len(phos_enz_codes)
#         if phos_enz_codes is not None
#         else 0,
#         "N_experiment": len(mi.experiments.experiments),
#         "N_flux_measurement": len(reaction_measurements),
#         "N_enzyme_measurement": len(enzyme_measurements),
#         "N_conc_measurement": len(mic_measurements),
#         "N_competitive_inhibitor": len(prior_dfs["inhibition_constant_priors"]),
#         "N_allosteric_inhibitor": len(prior_dfs["tense_dissociation_constant_priors"]),
#         "N_allosteric_activator": len(
#             prior_dfs["relaxed_dissociation_constant_priors"]
#         ),
#         "N_allosteric_enzyme": len(prior_dfs["transfer_constant_priors"]),
#         "N_drain": len(drain_codes) if drain_codes is not None else 0,
#         "unbalanced_mic_ix": list(unbalanced_mic_codes.values()),
#         "balanced_mic_ix": list(balanced_mic_codes.values()),
#         "experiment_yconc": (
#             mic_measurements["experiment_id"].map(experiment_codes).values
#         ),
#         "mic_ix_yconc": mic_measurements["target_id"].map(mic_codes).values,
#         "yconc": mic_measurements["value"].values,
#         "sigma_conc": mic_measurements["uncertainty"].values,
#         "experiment_yflux": (
#             reaction_measurements["experiment_id"].map(experiment_codes).values
#         ),
#         "reaction_yflux": (
#             reaction_measurements["target_id"].map(reaction_codes).values
#         ),
#         "yflux": reaction_measurements["value"].values,
#         "sigma_flux": reaction_measurements["uncertainty"].values,
#         "experiment_yenz": (
#             enzyme_measurements["experiment_id"].map(experiment_codes).values
#         ),
#         "enzyme_yenz": (enzyme_measurements["target_id"].map(enzyme_codes).values),
#         "yenz": enzyme_measurements["value"].values,
#         "sigma_enz": enzyme_measurements["uncertainty"].values,
#         "prior_loc_formation_energy": prior_dfs["formation_energy_priors"][
#             "location"
#         ].values,
#         "prior_scale_formation_energy": prior_dfs["formation_energy_priors"][
#             "scale"
#         ].values,
#         "prior_loc_kcat": prior_dfs["kcat_priors"]["location"].values,
#         "prior_scale_kcat": prior_dfs["kcat_priors"]["scale"].values,
#         "prior_loc_phos_kcat": prior_dfs["phos_kcat_priors"]["location"].values,
#         "prior_scale_phos_kcat": prior_dfs["phos_kcat_priors"]["scale"].values,
#         "prior_loc_km": prior_dfs["km_priors"]["location"].values,
#         "prior_scale_km": prior_dfs["km_priors"]["scale"].values,
#         "prior_loc_ki": prior_dfs["inhibition_constant_priors"]["location"].values,
#         "prior_scale_ki": prior_dfs["inhibition_constant_priors"]["scale"].values,
#         "prior_loc_diss_t": prior_dfs["tense_dissociation_constant_priors"][
#             "location"
#         ].values,
#         "prior_scale_diss_t": prior_dfs["tense_dissociation_constant_priors"][
#             "scale"
#         ].values,
#         "prior_loc_diss_r": prior_dfs["relaxed_dissociation_constant_priors"][
#             "location"
#         ].values,
#         "prior_scale_diss_r": prior_dfs["relaxed_dissociation_constant_priors"][
#             "scale"
#         ].values,
#         "prior_loc_tc": prior_dfs["transfer_constant_priors"]["location"].values,
#         "prior_scale_tc": prior_dfs["transfer_constant_priors"]["scale"].values,
#         "prior_loc_unbalanced": prior_loc_unb.values,
#         "prior_scale_unbalanced": prior_scale_unb.values,
#         "prior_loc_enzyme": prior_loc_enzyme.values,
#         "prior_scale_enzyme": prior_scale_enzyme.values,
#         "prior_loc_phos_conc": prior_loc_phos_conc.values,
#         "prior_scale_phos_conc": prior_scale_phos_conc.values,
#         "prior_loc_drain": prior_loc_drain.values.tolist(),
#         "prior_scale_drain": prior_scale_drain.values.tolist(),
#         "S_enz": enzyme_stoic.T.values,
#         "S_drain": drain_stoic.T.values.tolist(),
#         "S_full": full_stoic.T.values.tolist(),
#         "S_phos_act": S_phos_act.T.values.tolist(),
#         "S_phos_inh": S_phos_inh.T.values.tolist(),
#         "water_stoichiometry": water_stoichiometry,
#         "mic_to_met": list(mic_to_met.values()),
#         "S_to_flux_map": stoic_map_to_flux.values,
#         "is_knockout": knockout_matrix.values,
#         "is_phos_knockout": phos_knockout_matrix.values,
#         "km_lookup": km_lookup.values,
#         "n_ci": n_modifier["ki"],
#         "n_ai": n_modifier["diss_t"],
#         "n_aa": n_modifier["diss_r"],
#         "ci_ix": prior_dfs["inhibition_constant_priors"]["mic_code"].values,
#         "ai_ix": prior_dfs["tense_dissociation_constant_priors"]["mic_code"].values,
#         "aa_ix": prior_dfs["relaxed_dissociation_constant_priors"]["mic_code"].values,
#         "subunits": subunits["subunits"].values,
#         "conc_init": conc_init.values,
#         "init_enzyme": enz_conc_init.values,
#         "rel_tol": ode_config["rel_tol"],
#         "abs_tol": ode_config["abs_tol"],
#         "max_num_steps": int(ode_config["max_num_steps"]),
#         "LIKELIHOOD": int(likelihood),
#         "timepoint": ode_config["timepoint"],
#     }

