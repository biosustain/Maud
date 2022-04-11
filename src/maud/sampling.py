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

import json
import os
import warnings
from typing import Optional

import cmdstanpy
import numpy as np
import pandas as pd
from cmdstanpy.model import CmdStanModel
from cmdstanpy.stanfit.mcmc import CmdStanMCMC
from cmdstanpy.stanfit.vb import CmdStanVB
from xarray.core.dataset import Dataset

from maud.analysis import load_infd, load_infd_fit
from maud.utils import codify, get_null_space, get_rref
from maud.data_model.maud_input import MaudInput


HERE = os.path.dirname(os.path.abspath(__file__))
INCLUDE_PATH = ""
DEFAULT_PRIOR_LOC_DRAIN = None
DEFAULT_PRIOR_SCALE_DRAIN = None
STAN_PROGRAM_RELATIVE_PATH = "model.stan"
PPC_PROGRAM_RELATIVE_PATH = "out_of_sample_model.stan"

DEFAULT_SAMPLE_CONFIG = {
    "iter_warmup": 5,
    "iter_sampling": 5,
    "chains": 2,
    "max_treedepth": 11,
    "show_progress": True,
    "step_size": 0.025,
    "adapt_delta": 0.99,
    "save_warmup": True,
    "threads_per_chain": 1,
}
DEFAULT_VARIATIONAL_CONFIG = {
    "algorithm": "meanfield",
    "output_samples": 10,
    "require_converged": True,
}
SIM_CONFIG = {
    "chains": 1,
    "fixed_param": True,
    "iter_warmup": 0,
    "show_progress": False,
    "threads_per_chain": 1,
}


def sample(mi: MaudInput, output_dir: str) -> CmdStanMCMC:
    """Sample from the posterior defined by mi.

    :param mi: a MaudInput object
    :param output_dir: a string specifying where to save the output.
    """
    model = get_cmdstan_model(mi.config.cpp_options, mi.config.stanc_options)
    set_up_output_dir(output_dir, mi)
    sample_args: dict = {
        "data": os.path.join(output_dir, "input_data.json"),
        "inits": os.path.join(output_dir, "inits.json"),
        "output_dir": output_dir
    }
    sample_args = {**sample_args, **DEFAULT_SAMPLE_CONFIG}
    if mi.config.cmdstanpy_config is not None:
        sample_args = {**sample_args, **mi.config.cmdstanpy_config}
    return model.sample(**sample_args)


def variational(mi: MaudInput, output_dir: str) -> CmdStanVB:
    """Do variational inference for the posterior defined by mi.

    :param mi: a MaudInput object
    :param output_dir: a string specifying where to save the output.
    """
    mi_options = (
        {} if mi.config.variational_options is None else mi.config.variational_options
    )
    model = get_cmdstan_model(mi.config.cpp_options, mi.config.stanc_options)
    set_up_output_dir(output_dir, mi)
    return model.variational(
        data=os.path.join(output_dir, "input_data.json"),
        inits=os.path.join(output_dir, "inits.json"),
        **{
            **DEFAULT_VARIATIONAL_CONFIG,
            **mi_options,
            **{"output_dir": output_dir},
        },
    )


def simulate(mi: MaudInput, output_dir: str, n: int) -> cmdstanpy.CmdStanMCMC:
    """Generate simulations from the prior mean.

    :param mi: a MaudInput object
    :param output_dir: a string specifying where to save the output.
    """
    model = get_cmdstan_model(mi.config.cpp_options, mi.config.stanc_options)
    set_up_output_dir(output_dir, mi)
    return model.sample(
        data=os.path.join(output_dir, "input_data.json"),
        inits=os.path.join(output_dir, "inits.json"),
        **{**SIM_CONFIG, **{"output_dir": output_dir, "iter_sampling": n}},
    )


def set_up_output_dir(output_dir: str, mi: MaudInput):
    """Write input data, inits and coords to the output directory."""
    input_filepath = os.path.join(output_dir, "input_data.json")
    inits_filepath = os.path.join(output_dir, "inits.json")
    coords_filepath = os.path.join(output_dir, "coords.json")
    input_data = get_input_data(mi)
    inits = {k: v.values for k, v in mi.inits.items()}
    cmdstanpy.utils.write_stan_json(input_filepath, input_data)
    cmdstanpy.utils.write_stan_json(inits_filepath, inits)
    with open(coords_filepath, "w") as f:
        json.dump(mi.stan_coords.__dict__, f)


def get_cmdstan_model(
    cpp_options: Optional[dict], stanc_options: Optional[dict]
) -> CmdStanModel:
    """Get a CmdStanModel object.

    :param cpp_options: a dictionary of c++ options. For a list of things you
    can set see <cmdstan-home>/cmdstan/stan/lib/stan_math/make/compiler_flags

    :param stanc_options: a dictionary of c++ options. For a list of things you
    can set see https://mc-stan.org/docs/2_29/stan-users-guide/stanc-args.html

    """

    if cpp_options is None:
        cpp_options = {}
    if stanc_options is None:
        stanc_options = {}
    stanc_options = {
        **{"include-paths": [os.path.join(HERE, INCLUDE_PATH)]},
        **stanc_options,
    }
    return cmdstanpy.CmdStanModel(
        stan_file=os.path.join(HERE, STAN_PROGRAM_RELATIVE_PATH),
        stanc_options=stanc_options,
        cpp_options=cpp_options,
    )


def _sample_given_config(
    mi: MaudInput, output_dir: str, config: dict
) -> cmdstanpy.CmdStanMCMC:
    """Call CmdStanModel.sample, having already specified all arguments.

    :param mi: a MaudInput object
    :param output_dir: a string specifying where to save the output.
    :param config: a dictionary of keyword arguments to CmdStanModel.sample.
    """

    input_filepath = os.path.join(output_dir, "input_data.json")
    inits_filepath = os.path.join(output_dir, "inits.json")
    coords_filepath = os.path.join(output_dir, "coords.json")
    input_data = get_input_data(mi)
    # inits = {k: v.values for k, v in mi.inits.items()}
    cmdstanpy.utils.write_stan_json(input_filepath, input_data)
    cmdstanpy.utils.write_stan_json(inits_filepath, mi.inits)
    with open(coords_filepath, "w") as f:
        json.dump(mi.stan_coords.__dict__, f)
    config["inits"] = inits_filepath
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


def generate_predictions(
    mi_oos: MaudInput,
    mi_train: MaudInput,
    csvs: List[str],
    output_dir: str,
):
    """Call CmdStanModel.sample for out of sample predictions.

    :param mi_oos: a MaudInput object defining the out-of-sample experiments
    :param mi_train: a MaudInput object defining the training experiments
    :param output_dir: a string specifying where to save the output.
    """
    config = {
        **DEFAULT_SAMPLE_CONFIG,
        **mi_oos.config.cmdstanpy_config,
        **{"output_dir": output_dir},
    }
    input_filepath = os.path.join(output_dir, "input_data.json")
    coords_filepath = os.path.join(output_dir, "coords.json")
    input_data = get_input_data(mi_oos)
    cmdstanpy.utils.write_stan_json(input_filepath, input_data)
    with open(coords_filepath, "w") as f:
        json.dump(mi_oos.stan_coords.__dict__, f)
    ppc_program_filepath = os.path.join(HERE, PPC_PROGRAM_RELATIVE_PATH)
    include_path = os.path.join(HERE, INCLUDE_PATH)
    cpp_options = {}
    stanc_options = {"include_paths": [include_path]}
    infd = load_infd(csvs, mi_train)
    kinetic_parameters = [
        "keq",
        "km",
        "kcat",
        "diss_t",
        "diss_r",
        "transfer_constant",
        "kcat_phos",
        "ki",
    ]
    if config["threads_per_chain"] != 1:
        cpp_options["STAN_THREADS"] = True
        os.environ["STAN_NUM_THREADS"] = str(config["threads_per_chain"])
    posterior = infd.get("posterior")
    assert isinstance(posterior, Dataset)
    chains = posterior.chain.to_series().values
    draws = posterior.draw.to_series().values
    gq_model = cmdstanpy.CmdStanModel(
        stan_file=ppc_program_filepath,
        stanc_options=stanc_options,
        cpp_options=cpp_options,
    )
    all_conc = []
    all_conc_enz = []
    all_flux = []
    for chain in chains:
        for draw in draws:
            inits = {
                par: posterior[par][chain][draw].to_series().values
                for par in kinetic_parameters
                if par in posterior.variables.keys()
            }
            gq_samples = gq_model.sample(
                inits=inits,
                iter_warmup=0,
                iter_sampling=1,
                data=input_filepath,
                fixed_param=True,
            )
            posterior_fit = load_infd_fit(gq_samples.runset.csv_files, mi_oos).get(
                "posterior"
            )
            assert isinstance(posterior_fit, Dataset)
            tmp_conc = posterior_fit.conc.to_dataframe().reset_index()
            tmp_conc["chain"] = chain
            tmp_conc["draw"] = draw
            tmp_conc_enz = posterior_fit.conc_enzyme.to_dataframe().reset_index()
            tmp_conc_enz["chain"] = chain
            tmp_conc_enz["draw"] = draw
            tmp_flux = posterior_fit.flux.to_dataframe().reset_index()
            tmp_flux["chain"] = chain
            tmp_flux["draw"] = draw
            all_conc.append(tmp_conc)
            all_conc_enz.append(tmp_conc_enz)
            all_flux.append(tmp_flux)
    draws = {}
    flux_df = pd.concat(all_flux)
    conc_df = pd.concat(all_conc)
    conc_enz_df = pd.concat(all_conc_enz)
    flux_df.to_csv(os.path.join(output_dir, "flux.csv"))
    conc_df.to_csv(os.path.join(output_dir, "conc.csv"))
    conc_enz_df.to_csv(os.path.join(output_dir, "conc_enzyme.csv"))


def validate_specified_fluxes(mi: MaudInput):
    """Check that appropriate fluxes have been measured.

    :param mi: A MaudInput object
    """
    S = get_stoichiometry(mi)
    balanced_mic_ix = mi.stan_coords.balanced_mics
    reactions = mi.stan_coords.reactions
    for exp in mi.stan_coords.experiments:
        measured_rxns = mi.measurements.yflux.reset_index()["target_id"].unique()
        flux_paths = get_null_space(S.loc[balanced_mic_ix].values)
        _, n_dof = np.shape(flux_paths)
        rref_flux_paths = np.matrix(get_rref(flux_paths.T))
        rref_flux_paths[np.abs(rref_flux_paths) < 1e-10] = 0
        flux_paths_df = pd.DataFrame(rref_flux_paths, columns=reactions)
        for _, flux_path in flux_paths_df.iterrows():
            if any(flux_paths[measured_rxns]) != 0:
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
