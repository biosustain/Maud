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
from typing import Dict

import cmdstanpy
import numpy as np
import pandas as pd

from maud import code_generation, io, utils
from maud.data_model import KineticModel, MaudInput


INCLUDE_PATH = "stan_code"
DEFAULT_PRIOR_LOC_UNBALANCED = 1
DEFAULT_PRIOR_SCALE_UNBALANCED = 4
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


def prior_check(
    data_path: str,
    n_samples: int,
    n_warmup: int,
    n_chains: int,
    n_cores: int,
    output_dir: str,
) -> cmdstanpy.CmdStanMCMC:
    """Sample from a prior distribution.

    :param data_path: A path to a toml file containing input data
    :param n_samples: Number of post-warmup samples
    :param n_warmup: Number of warmup samples
    :param n_chains: Number of MCMC chains to run
    :param n_cores: Number of cores to try and use
    :param: output_dir: Directory to save output
    """
    model_name = os.path.splitext(os.path.basename(data_path))[0]
    input_filepath = os.path.join(output_dir, f"input_data_prior_{model_name}.json")

    here = os.path.dirname(os.path.abspath(__file__))

    mi = io.load_maud_input_from_toml(data_path)

    input_data = get_input_data(mi)

    cmdstanpy.utils.jsondump(input_filepath, input_data)

    stan_program_filename = f"prior_model_{model_name}.stan"
    stan_program_filepath = os.path.join(output_dir, stan_program_filename)
    exe_file_path = stan_program_filepath[:-5]
    stan_code = code_generation.create_stan_program(mi, "prior_check")
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

    model_fit = model.sample(
        data=input_filepath,
        chains=n_chains,
        cores=4,
        iter_sampling=n_samples,
        output_dir=output_dir,
        iter_warmup=n_warmup,
        max_treedepth=15,
        save_warmup=True,
        show_progress=True,
    )

    reaction_measurements = pd.DataFrame(
            [
                [exp.id, meas.target_id, meas.value, meas.uncertainty]
                for exp in mi.experiments.values()
                for meas in exp.measurements['reaction'].values()
            ],
            columns=["experiment_id", "target_id", "value", "uncertainty"],
        )

    for i, row in enumerate(reaction_measurements.iterrows()):
        if row[1]['value'] < model_fit.summary().loc[f'yflux_sim[{i+1}]', '5%']:
            print(f"""{row[1]['target_id']} in experiment {row[1]['experiment_id']}
                had a flux specified that was smaller than the prior flux\n
                flux = {row[1]['value']} \t\t limit = {model_fit.summary().loc[f'yflux_sim[{i+1}]', '5%']}""")

        if row[1]['value'] > model_fit.summary().loc[f'yflux_sim[{i+1}]', '95%']:
            print(f"""{row[1]['target_id']} in experiment {row[1]['experiment_id']}
                had a flux specified that was larger than the prior flux\n
                flux = {row[1]['value']} \t\t limit = {model_fit.summary().loc[f'yflux_sim[{i+1}]', '95%']}""")

    return


def get_input_data(
    mi: MaudInput,
) -> dict:
    """Put a MaudInput and some config numbers into a Stan-friendly dictionary.

    :param mi: a MaudInput object
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
    mic_shape = len(mi.experiments), len(mic_codes)
    prior_loc_conc = np.full(mic_shape, DEFAULT_PRIOR_LOC_UNBALANCED)
    prior_scale_conc = np.full(mic_shape, DEFAULT_PRIOR_SCALE_UNBALANCED)
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
        "prior_loc_conc": prior_loc_conc,
        "prior_scale_conc": prior_scale_conc,
        "prior_loc_enzyme": prior_loc_enzyme,
        "prior_scale_enzyme": prior_scale_enzyme,
    }
