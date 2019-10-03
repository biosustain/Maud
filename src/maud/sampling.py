# Copyright (c) 2019, Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Code for sampling from a posterior distribution."""

import os

import arviz
import cmdstanpy
import numpy as np
import pandas as pd

from maud import code_generation, io, utils
from scipy.linalg import null_space as null_space


RELATIVE_PATHS = {
    "stan_includes": "stan_code",
    "stan_autogen": "stan_code/autogen",
    "stan_records": "../../data/stan_records",
    "data_out": "../../data/out",
}


def get_full_stoichiometry(kinetic_model, enzyme_codes, metabolite_codes):
    """Gets the full stoichiometric matrix for each isoenzyme

    :param kinetic_model: A Kinetic Model object
    :param enzyme_codes: the codified enzyme codes
    """
    S = pd.DataFrame(index=enzyme_codes, columns=metabolite_codes)
    for _, rxn in kinetic_model.reactions.items():
        for enz_id, _ in rxn.enzymes.items():
            for met, stoic in rxn.stoichiometry.items():
                S.loc[enz_id, met] = stoic
    S.fillna(0, inplace=True)
    return S

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
):
    """Get samples from a posterior distribution.

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

    # define input data
    mi = io.load_maud_input_from_toml(data_path)
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
    # stan codes
    experiment_codes = utils.codify(mi.experiments.keys())
    reaction_codes = utils.codify(reactions.keys())
    enzyme_codes = utils.codify(enzymes.keys())
    metabolite_codes = utils.codify(metabolites.keys())
    full_stoic = get_full_stoichiometry(mi.kinetic_model, enzyme_codes, metabolite_codes)
    flux_nullspace = null_space(np.transpose(np.matrix(full_stoic)))
    flux_nullspace[(flux_nullspace < 10e-10) & (flux_nullspace > -10e-10)] = 0

    if not flux_nullspace.any():
        wegschneider_mat = np.identity(len(enzyme_codes))
        n_loops = len(enzyme_codes)
    else:
        wegschneider_mat = null_space(np.transpose(flux_nullspace))
        n_loops = np.shape(wegschneider_mat)[1]
    
    input_data = {
        "N_balanced": len(balanced_metabolites),
        "N_unbalanced": len(unbalanced_metabolites),
        "N_kinetic_parameter": len(kinetic_parameter_priors),
        "N_reaction": len(reactions),
        "N_enzyme": len(enzymes),
        "N_experiment": len(mi.experiments),
        "N_flux_measurement": len(reaction_measurements),
        "N_conc_measurement": len(metabolite_measurements),
        "N_loops": n_loops,
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
        "ln_equilibrium_basis": wegschneider_mat,
        "rel_tol": rel_tol,
        "f_tol": f_tol,
        "max_steps": max_steps,
        "LIKELIHOOD": likelihood,
    }

    # dump input data
    input_file = os.path.join(paths["stan_records"], f"input_data_{model_name}.json")
    cmdstanpy.utils.jsondump(input_file, input_data)

    # compile model if necessary
    stan_code = code_generation.create_stan_program(mi, "inference", time_step)
    stan_file = os.path.join(
        paths["stan_autogen"], f"inference_model_{model_name}.stan"
    )
    exe_file = stan_file[:-5]
    no_need_to_compile = os.path.exists(exe_file) and utils.match_string_to_file(
        stan_code, stan_file
    )
    if no_need_to_compile:
        model = cmdstanpy.Model(stan_file=stan_file, exe_file=exe_file)
    else:
        with open(stan_file, "w") as f:
            f.write(stan_code)
        model = cmdstanpy.Model(stan_file)
        model.compile(include_paths=[paths["stan_includes"]], overwrite=True)

    # draw samples
    csv_output_file = os.path.join(paths["data_out"], f"output_{model_name}.csv")

    fit = model.sample(
        data=input_file,
        cores=4,
        chains=n_chains,
        csv_basename=csv_output_file,
        sampling_iters=n_samples,
        warmup_iters=n_warmup,
        max_treedepth=15,
        adapt_delta=0.8,
        save_warmup=True,
        inits={
            "kinetic_parameter": np.exp(kinetic_parameter_priors["location"]).T.values,
            "unbalanced": np.exp(prior_loc_unb).T.values,
            "enzyme_concentration": np.exp(prior_loc_enzyme).T.values,
        },
    )

    infd_posterior = arviz.from_cmdstanpy(
        posterior=fit,
        posterior_predictive=["yflux_sim", "yconc_sim"],
        observed_data={
            "yflux_sim": input_data["yflux"],
            "yconc_sim": input_data["yconc"],
        },
        coords={
            "reactions": list(reaction_codes.keys()),
            "metabolites": list(metabolite_codes.keys()),
            "experiments": list(experiment_codes.keys()),
            "enzymes": list(enzyme_codes.keys()),
            "kinetic_parameter_names": kinetic_parameter_priors["id"].tolist(),
            "reaction_measurements": [
                str(experiment_codes[row["experiment_id"]])
                + "_"
                + str(reaction_codes[row["target_id"]])
                for _, row in reaction_measurements.iterrows()
            ],
            "metabolite_measurements": [
                str(experiment_codes[row["experiment_id"]])
                + "_"
                + str(metabolite_codes[row["target_id"]])
                for _, row in metabolite_measurements.iterrows()
            ],
        },
        dims={
            "conc": ["experiments", "metabolites"],
            "flux": ["experiments", "reactions"],
            "kinetic_parameter": ["kinetic_parameter_names"],
            "enzyme_concentration": ["experiments", "enzymes"],
            "yconc_sim": ["metabolite_measurements"],
            "yflux_sim": ["reaction_measurements"],
        },
    )

    infd_posterior.to_netcdf(
        os.path.join(paths["data_out"], f"model_inference_{model_name}.nc")
    )

    return fit
