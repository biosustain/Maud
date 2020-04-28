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
from copy import deepcopy
from typing import Dict

import cmdstanpy
import numpy as np
import pandas as pd

from maud import code_generation, io, utils
from maud.data_model import KineticModel, MaudInput


INCLUDE_PATH = "stan_code"
DEFAULT_PRIOR_LOC_UNBALANCED = 0.1
DEFAULT_PRIOR_SCALE_UNBALANCED = 3
DEFAULT_PRIOR_LOC_ENZYME = 0.1
DEFAULT_PRIOR_SCALE_ENZYME = 3


def get_full_stoichiometry(
    kinetic_model: KineticModel,
    enzyme_codes: Dict[str, int],
    metabolite_codes: Dict[str, int],
):
    """Get full stoichiometric matrix for each isoenzyme.

    :param kinetic_model: A Kinetic Model object
    :param enzyme_codes: the codified enzyme codes
    :param metabolite_codes: the codified metabolite codes
    """
    S = pd.DataFrame(index=enzyme_codes, columns=metabolite_codes)
    for _, rxn in kinetic_model.reactions.items():
        for enz_id, _ in rxn.enzymes.items():
            for met, stoic in rxn.stoichiometry.items():
                S.loc[enz_id, met] = stoic
    S.fillna(0, inplace=True)
    return S


def get_knockout_matrix(mi: MaudInput):
    """Get binary experiment, enzyme matrix, 1 if enyzme present, 0 if not.

    :param mi: a MaudInput object
    """
    experiment_codes = mi.stan_codes["experiment"]
    enzyme_codes = mi.stan_codes["enzyme"]

    enzyme_knockout_matrix = pd.DataFrame(
        1, index=np.arange(len(experiment_codes)), columns=np.arange(len(enzyme_codes))
    )

    for exp_name, exp in mi.experiments.items():
        if exp.knockouts:
            for enz in exp.knockouts:
                enzyme_knockout_matrix.loc[
                    experiment_codes[exp_name] - 1, enzyme_codes[enz] - 1
                ] = 0

    return enzyme_knockout_matrix


def sample(
    data_path: str,
    f_tol: float,
    rel_tol: float,
    max_steps: int,
    likelihood: int,
    n_samples: int,
    n_warmup: int,
    n_chains: int,
    n_cores: int,
    timepoint: float,
    output_dir: str,
) -> cmdstanpy.CmdStanMCMC:
    """Sample from a posterior distribution.

    :param data_path: A path to a toml file containing input data
    :param f_tol: Sets algebra solver's f_tol control parameter
    :param rel_tol: Sets algebra solver's rel_tol control parameter
    :param max_steps: Sets algebra solver's max_steps control parameter
    :param likelihood: Set to zero if you don't want the model to include information
    from experimental data.
    :param n_samples: Number of post-warmup samples
    :param n_warmup: Number of warmup samples
    :param n_chains: Number of MCMC chains to run
    :param n_cores: Number of cores to try and use
    :param time_step: Amount of time for the ode solver to simulate in order to compare
    initial state with evolved state
    :param: output_dir: Directory to save output
    """
    model_name = os.path.splitext(os.path.basename(data_path))[0]
    input_filepath = os.path.join(output_dir, f"input_data_{model_name}.json")

    here = os.path.dirname(os.path.abspath(__file__))

    mi = io.load_maud_input_from_toml(data_path)

    input_data = get_input_data(mi, f_tol, rel_tol, max_steps, likelihood, timepoint)
    init_cond = get_initial_conditions(input_data, mi)

    cmdstanpy.utils.jsondump(input_filepath, input_data)

    stan_program_filename = f"inference_model_{model_name}.stan"
    stan_program_filepath = os.path.join(output_dir, stan_program_filename)
    exe_file_path = stan_program_filepath[:-5]
    stan_code = code_generation.create_stan_program(mi, "inference")
    exe_file_exists = os.path.exists(exe_file_path)
    change_in_stan_code = not utils.match_string_to_file(
        stan_code, stan_program_filepath
    )
    need_to_overwrite = (not exe_file_exists) or change_in_stan_code
    if need_to_overwrite:
        with open(stan_program_filepath, "w") as f:
            f.write(stan_code)
        for p in [exe_file_path, exe_file_path + ".o", exe_file_path + ".hpp"]:
            if os.path.exists(p):
                os.remove(p)
    include_path = os.path.join(here, INCLUDE_PATH)
    model = cmdstanpy.CmdStanModel(
        stan_file=stan_program_filepath, stanc_options={"include_paths": [include_path]}
    )
    return model.sample(
        data=input_filepath,
        chains=n_chains,
        cores=4,
        iter_sampling=n_samples,
        output_dir=output_dir,
        iter_warmup=n_warmup,
        max_treedepth=15,
        save_warmup=True,
        inits=init_cond,
        show_progress=True,
    )


def get_input_data(
    mi: MaudInput,
    f_tol: float,
    rel_tol: float,
    max_steps: int,
    likelihood: int,
    timepoint: float,
) -> dict:
    """Put a MaudInput and some config numbers into a Stan-friendly dictionary.

    :param mi: a MaudInput object
    :param f_tol: Sets algebra solver's f_tol control parameter
    :param rel_tol: Sets algebra solver's rel_tol control parameter
    :param max_steps: Sets algebra solver's max_steps control parameter
    :param likelihood: Set to zero if you don't want the model to include information
    """
    metabolites = mi.kinetic_model.metabolites
    mics = mi.kinetic_model.mics
    reactions = mi.kinetic_model.reactions
    reaction_codes = mi.stan_codes["reaction"]
    enzyme_codes = mi.stan_codes["enzyme"]
    kp_codes = mi.stan_codes["kinetic_parameter"]
    experiment_codes = mi.stan_codes["experiment"]
    met_codes = mi.stan_codes["metabolite"]
    mic_codes = mi.stan_codes["metabolite_in_compartment"]
    balanced_mic_codes = mi.stan_codes["balanced_mic"]
    unbalanced_mic_codes = mi.stan_codes["unbalanced_mic"]
    mic_to_met = {
        mic_codes[mic.id]: met_codes[mic.metabolite_id] for mic in mics.values()
    }
    enzymes = {k: v for r in reactions.values() for k, v in r.enzymes.items()}
    full_stoic = get_full_stoichiometry(mi.kinetic_model, enzyme_codes, mic_codes)
    # priors
    prior_loc_kp = [mi.priors[k].location for k in kp_codes.keys()]
    prior_scale_kp = [mi.priors[k].scale for k in kp_codes.keys()]
    unb_shape = len(mi.experiments), len(unbalanced_mic_codes)
    prior_loc_unb = np.full(unb_shape, DEFAULT_PRIOR_LOC_UNBALANCED)
    prior_scale_unb = np.full(unb_shape, DEFAULT_PRIOR_SCALE_UNBALANCED)
    for p in mi.priors.values():
        if p.target_type == "metabolite_concentration":
            ix = [experiment_codes[mic_codes[p.target_id - 1], p.experiment_id] - 1]
            prior_loc_unb[ix] = p.location
            prior_scale_unb[ix] = p.scale
    prior_loc_formation_energy = [
        mi.priors[k + "_formation_energy"].location for k in met_codes.keys()
    ]
    prior_scale_formation_energy = [
        mi.priors[k + "_formation_energy"].scale for k in met_codes.keys()
    ]
    enzyme_shape = len(mi.experiments), len(enzymes)
    prior_loc_enzyme = np.full(enzyme_shape, DEFAULT_PRIOR_LOC_ENZYME)
    prior_scale_enzyme = np.full(enzyme_shape, DEFAULT_PRIOR_SCALE_ENZYME)
    # measurements
    mic_measurements, reaction_measurements, enzyme_measurements = (
        pd.DataFrame(
            [
                [exp.id, meas.target_id, meas.value, meas.uncertainty]
                for exp in mi.experiments.values()
                for meas in exp.measurements[measurement_type].values()
            ],
            columns=["experiment_id", "target_id", "value", "uncertainty"],
        )
        for measurement_type in ["metabolite", "reaction", "enzyme"]
    )

    balanced_init = pd.DataFrame(
        0.01, index=experiment_codes.values(), columns=balanced_mic_codes.values()
    )
    for _i, row in mic_measurements.iterrows():
        if row["target_id"] in balanced_mic_codes.keys():
            row_ix = experiment_codes[row["experiment_id"]]
            column_ix = balanced_mic_codes[row["target_id"]]
            balanced_init.loc[row_ix, column_ix] = row["value"]

    knockout_matrix = get_knockout_matrix(mi=mi)

    return {
        "N_mic": len(mics),
        "N_unbalanced": len(unbalanced_mic_codes),
        "N_kinetic_parameters": len(kp_codes),
        "N_reaction": len(reactions),
        "N_enzyme": len(enzymes),
        "N_experiment": len(mi.experiments),
        "N_flux_measurement": len(reaction_measurements),
        "N_enzyme_measurement": len(enzyme_measurements),
        "N_conc_measurement": len(mic_measurements),
        "N_metabolite": len(metabolites),
        "stoichiometric_matrix": full_stoic.T.values,
        "metabolite_ix_stoichiometric_matrix": list(mic_to_met.values()),
        "experiment_yconc": (
            mic_measurements["experiment_id"].map(experiment_codes).values
        ),
        "mic_ix_yconc": mic_measurements["target_id"].map(mic_codes).values,
        "balanced_mic_ix": list(balanced_mic_codes.values()),
        "unbalanced_mic_ix": list(unbalanced_mic_codes.values()),
        "yconc": mic_measurements["value"].values,
        "sigma_conc": mic_measurements["uncertainty"].values,
        "experiment_yflux": (
            reaction_measurements["experiment_id"].map(experiment_codes).values
        ),
        "reaction_yflux": (
            reaction_measurements["target_id"].map(reaction_codes).values
        ),
        "yflux": reaction_measurements["value"].values,
        "sigma_flux": reaction_measurements["uncertainty"].values,
        "experiment_yenz": (
            enzyme_measurements["experiment_id"].map(experiment_codes).values
        ),
        "enzyme_yenz": (enzyme_measurements["target_id"].map(enzyme_codes).values),
        "yenz": enzyme_measurements["value"].values,
        "sigma_enz": enzyme_measurements["uncertainty"].values,
        "prior_loc_formation_energy": prior_loc_formation_energy,
        "prior_scale_formation_energy": prior_scale_formation_energy,
        "prior_loc_kinetic_parameters": prior_loc_kp,
        "prior_scale_kinetic_parameters": prior_scale_kp,
        "prior_loc_unbalanced": prior_loc_unb,
        "prior_scale_unbalanced": prior_scale_unb,
        "prior_loc_enzyme": prior_loc_enzyme,
        "prior_scale_enzyme": prior_scale_enzyme,
        "conc_init": balanced_init.values,
        "knockout_enzymes": knockout_matrix.values,
        "rtol": rel_tol,
        "ftol": f_tol,
        "steps": max_steps,
        "LIKELIHOOD": likelihood,
        "timepoint": timepoint,
    }


def validate_init_enzyme(mi):
    """Initilize enyme concentrations above their measured flux.

    :param mi: a MaudInput object
    """
    enzyme_validation = pd.DataFrame(
        columns=[
            "reaction",
            "enzyme",
            "experiment",
            "enzyme_concentration",
            "kcat_prior",
            "reaction_flux",
            "scaling_factor",
            "init_enzyme",
        ]
    )

    rxn_enz_dict = {
        rxn_name: list(rxn_data.enzymes.keys())
        for rxn_name, rxn_data in mi.kinetic_model.reactions.items()
    }

    count = 0

    for exp_name, exp in mi.experiments.items():
        for rxn, flux in exp.measurements["reaction"].items():
            tmp_enz_kcat_array = []
            Vmax = 0
            for enz in rxn_enz_dict[rxn]:
                try:
                    tmp_enz_kcat_array.append(
                        [
                            enz,
                            exp.measurements["enzyme"][enz].value,
                            mi.priors[f"{enz}_Kcat1"].location,
                        ]
                    )
                except KeyError:
                    print(
                        f"The enzyme: {enz} wasn't measured in experiment: {exp_name}"
                    )
                    tmp_enz_kcat_array.append(
                        [
                            enz,
                            DEFAULT_PRIOR_LOC_ENZYME,
                            mi.priors[f"{enz}_Kcat1"].location,
                        ]
                    )
                Vmax += tmp_enz_kcat_array[-1][1] * tmp_enz_kcat_array[-1][2]
            if Vmax < np.abs(flux.value):
                scaling_factor = 1.5 * flux.value / Vmax
            else:
                scaling_factor = 1

            for enz in tmp_enz_kcat_array:
                enzyme_validation.loc[count, "reaction"] = rxn
                enzyme_validation.loc[count, "enzyme"] = enz[0]
                enzyme_validation.loc[count, "experiment"] = exp_name
                enzyme_validation.loc[count, "enzyme_concentration"] = enz[1]
                enzyme_validation.loc[count, "kcat_prior"] = enz[2]
                enzyme_validation.loc[count, "reaction_flux"] = np.abs(flux.value)
                enzyme_validation.loc[count, "scaling_factor"] = scaling_factor
                enzyme_validation.loc[count, "init_enzyme"] = np.abs(
                    scaling_factor * enz[1]
                )
                count += 1
    return enzyme_validation


def get_initial_conditions(input_data, mi):
    """Specify parameters' initial conditions."""
    experiment_codes = mi.stan_codes["experiment"]
    enzyme_codes = mi.stan_codes["enzyme"]
    init_conc_unb = deepcopy(input_data["prior_loc_unbalanced"])

    init_unbalanced = pd.DataFrame(
        init_conc_unb,
        index=range(1, input_data["N_experiment"] + 1),
        columns=input_data["unbalanced_mic_ix"],
    )
    for exp_ix, mic_ix, measurement in zip(
        input_data["experiment_yconc"], input_data["mic_ix_yconc"], input_data["yconc"]
    ):
        if mic_ix in input_data["unbalanced_mic_ix"]:
            init_unbalanced.loc[exp_ix, mic_ix] = measurement

    init_enzyme_df = validate_init_enzyme(mi)

    init_enzyme = pd.DataFrame(
        index=range(1, input_data["N_experiment"] + 1),
        columns=sorted(pd.unique(input_data["enzyme_yenz"])),
    )

    for _, row in init_enzyme_df.iterrows():
        init_enzyme.loc[
            experiment_codes[row["experiment"]], enzyme_codes[row["enzyme"]]
        ] = row["init_enzyme"]

    init_enzyme.fillna(DEFAULT_PRIOR_LOC_ENZYME, inplace=True)

    return {
        "kinetic_parameters": input_data["prior_loc_kinetic_parameters"],
        "conc_unbalanced": init_unbalanced.values,
        "enzyme_concentration": init_enzyme.values,
        "formation_energy": input_data["prior_loc_formation_energy"],
    }
