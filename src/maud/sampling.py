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

import cmdstanpy
import pandas as pd

from maud import code_generation, io, utils
from maud.data_model import MaudInput


RELATIVE_PATHS = {
    "stan_includes": "stan_code",
    "autogen": "stan_code/autogen",
    "stan_records": "../../data/stan_records",
    "data_out": "../../data/out",
}


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
    time_step: float,
) -> cmdstanpy.StanFit:
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
    """

    model_name = os.path.splitext(os.path.basename(data_path))[0]
    here = os.path.dirname(os.path.abspath(__file__))
    paths = {k: os.path.join(here, v) for k, v in RELATIVE_PATHS.items()}

    mi = io.load_maud_input_from_toml(data_path)

    input_data = get_input_data(mi, f_tol, rel_tol, max_steps, likelihood)
    input_file = os.path.join(paths["stan_records"], f"input_data_{model_name}.json")
    cmdstanpy.utils.jsondump(input_file, input_data)

    stan_file = os.path.join(paths["autogen"], f"inference_model_{model_name}.stan")
    stan_code = code_generation.create_stan_program(mi, "inference", time_step)
    no_exe_file = not os.path.exists(stan_file[:-5])
    change_in_stan_code = not utils.match_string_to_file(stan_code, stan_file)
    need_to_overwrite = no_exe_file or change_in_stan_code
    model = cmdstanpy.Model(stan_file=stan_file)
    model.compile(include_paths=[paths["stan_includes"]], overwrite=need_to_overwrite)

    return model.sample(
        data=input_file,
        cores=4,
        chains=n_chains,
        csv_basename=os.path.join(paths["data_out"], f"output_{model_name}.csv"),
        sampling_iters=n_samples,
        warmup_iters=n_warmup,
        max_treedepth=15,
        save_warmup=True,
    )


def get_input_data(
    mi: MaudInput, f_tol: float, rel_tol: float, max_steps: int, likelihood: int
) -> dict:
    """Put a MaudInput and some config numbers into a Stan-friendly dictionary.

    :param mi: a MaudInput object
    :param f_tol: Sets algebra solver's f_tol control parameter
    :param rel_tol: Sets algebra solver's rel_tol control parameter
    :param max_steps: Sets algebra solver's max_steps control parameter
    :param likelihood: Set to zero if you don't want the model to include information
    """
    prior_df = pd.DataFrame.from_records(
        [
            [p.id, p.experiment_id, p.target_id, p.location, p.scale, p.target_type]
            for p in mi.priors.values()
        ],
        columns=[
            "id",
            "experiment_id",
            "target_id",
            "location",
            "scale",
            "target_type",
        ],
    )
    metabolites = mi.kinetic_model.metabolites
    reactions = mi.kinetic_model.reactions
    enzymes = {k: v for r in reactions.values() for k, v in r.enzymes.items()}
    balanced_metabolites = {k: v for k, v in metabolites.items() if v.balanced}
    unbalanced_metabolites = {k: v for k, v in metabolites.items() if not v.balanced}
    unbalanced_metabolite_priors, kinetic_parameter_priors, enzyme_priors = (
        prior_df.loc[lambda df: df["target_type"] == target_type]
        for target_type in ["unbalanced_metabolite", "kinetic_parameter", "enzyme"]
    )
    prior_loc_unb, prior_loc_enzyme, prior_scale_unb, prior_scale_enzyme = (
        df.set_index(["target_id", "experiment_id"])[col].unstack()
        for col in ["location", "scale"]
        for df in [unbalanced_metabolite_priors, enzyme_priors]
    )
    metabolite_measurements, reaction_measurements = (
        pd.DataFrame(
            [
                [exp.id, meas.target_id, meas.value, meas.uncertainty]
                for exp in mi.experiments.values()
                for meas in exp.measurements[measurement_type].values()
            ],
            columns=["experiment_id", "target_id", "value", "uncertainty"],
        )
        for measurement_type in ["metabolite", "reaction"]
    )
    experiment_codes = utils.codify(mi.experiments.keys())
    reaction_codes = utils.codify(reactions.keys())
    metabolite_codes = utils.codify(metabolites.keys())
    return {
        "N_balanced": len(balanced_metabolites),
        "N_unbalanced": len(unbalanced_metabolites),
        "N_kinetic_parameter": len(kinetic_parameter_priors),
        "N_reaction": len(reactions),
        "N_enzyme": len(enzymes),
        "N_experiment": len(mi.experiments),
        "N_flux_measurement": len(reaction_measurements),
        "N_conc_measurement": len(metabolite_measurements),
        "experiment_yconc": (
            metabolite_measurements["experiment_id"].map(experiment_codes).values
        ),
        "metabolite_yconc": (
            metabolite_measurements["target_id"].map(metabolite_codes).values
        ),
        "yconc": metabolite_measurements["value"].values,
        "sigma_conc": metabolite_measurements["uncertainty"].values,
        "experiment_yflux": (
            reaction_measurements["experiment_id"].map(experiment_codes).values
        ),
        "reaction_yflux": (
            reaction_measurements["target_id"].map(reaction_codes).values
        ),
        "yflux": reaction_measurements["value"].values,
        "sigma_flux": reaction_measurements["uncertainty"].values,
        "prior_loc_kinetic_parameter": kinetic_parameter_priors["location"].values,
        "prior_scale_kinetic_parameter": kinetic_parameter_priors["scale"].values,
        "prior_loc_unbalanced": prior_loc_unb.values,
        "prior_scale_unbalanced": prior_scale_unb.values,
        "prior_loc_enzyme": prior_loc_enzyme.values,
        "prior_scale_enzyme": prior_scale_enzyme.values,
        "balanced_guess": [1.0 for m in range(len(balanced_metabolites))],
        "rel_tol": rel_tol,
        "f_tol": f_tol,
        "max_steps": max_steps,
        "LIKELIHOOD": likelihood,
    }
