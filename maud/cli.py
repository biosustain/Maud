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

"""Functions that are exposed to the command line interface live here."""

import json
import os
import shutil
from datetime import datetime

import arviz as az
import click
import importlib_resources
from stanio import write_stan_json

from maud.data.example_inputs import linear, methionine
from maud.getting_idatas import get_idata
from maud.loading_maud_inputs import load_maud_input
from maud.running_stan import (
    laplace,
    optimize,
    pathfinder,
    predict,
    sample,
    simulate,
    variational,
)

AVAILABLE_EXAMPLE_INPUTS = {"linear": linear, "methionine": methionine}


@click.group()
@click.help_option("--help", "-h")
def cli():
    """Use Maud's command line interface."""
    pass


@cli.command("load-input")
@click.option("--output_dir", default=".", help="Where to save the folder.")
@click.argument("input_name", type=click.STRING)
def load_input_command(input_name, output_dir):
    """Extract one of Maud's example inputs.

    Available inputs: "linear", "methionine"
    """
    if input_name not in AVAILABLE_EXAMPLE_INPUTS.keys():
        raise ValueError(
            f"{input_name} not available. Avaliable inputs:\n\t"
            + "\n\t".join(AVAILABLE_EXAMPLE_INPUTS.keys())
        )
    data = AVAILABLE_EXAMPLE_INPUTS[input_name]
    path = importlib_resources.files(data)._paths[0]
    shutil.copytree(path, os.path.join(output_dir, input_name))


@cli.command("sample")
@click.option("--output_dir", default=".", help="Where to save Maud's output")
@click.argument(
    "data_path",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
)
def sample_command(data_path, output_dir):
    """Generate MCMC samples given a user input directory.

    This function creates a new directory in output_dir with a name starting
    with "maud_output". It first copies the directory at data_path into the new
    this directory at new_dir/user_input, then runs the running_stan.sample
    function to write samples in new_dir/samples. Finally it prints the results
    of cmdstanpy's diagnose and summary methods.


    Run the sample function as a click command.
    """
    click.echo(do_sample(data_path, output_dir))


def do_sample(data_path, output_dir):
    """Generate MCMC samples given a user input directory.

    This function creates a new directory in output_dir with a name starting
    with "maud_output". It first copies the directory at data_path into the new
    this directory at new_dir/user_input, then runs the running_stan.sample
    function to write samples in new_dir/samples. Finally it prints the results
    of cmdstanpy's diagnose and summary methods.

    """
    mi = load_maud_input(data_path)
    now = datetime.now().strftime("%Y%m%d%H%M%S")
    output_name = f"maud_output-{mi.config.name}-{now}"
    output_path = os.path.join(output_dir, output_name)
    samples_path = os.path.join(output_path, "samples")
    ui_dir = os.path.join(output_path, "user_input")
    print("Creating output directory: " + output_path)
    os.mkdir(output_path)
    os.mkdir(samples_path)
    print(f"Copying user input from {data_path} to {ui_dir}")
    shutil.copytree(data_path, ui_dir)
    stanfit = sample(mi, samples_path)
    print(stanfit.diagnose())
    print(stanfit.summary())
    if stanfit.metric is not None:
        for i, chain_metric in enumerate(stanfit.metric):
            chain_metric_dict = {"inv_metric": chain_metric.tolist()}
            with open(
                os.path.join(samples_path, f"metric-{str(i)}.json"), "w"
            ) as f:
                json.dump(chain_metric_dict, f)
    idata = get_idata(stanfit.runset.csv_files, mi, "train")
    idata.to_json(os.path.join(output_path, "idata.json"))
    return output_path


@cli.command("predict")
@click.argument(
    "data_path",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
)
def predict_command(data_path):
    """Generate MCMC samples given a Maud output folder at train_path.

    This function creates a new directory in output_dir with a name starting
    with "maud-predict-output". It first copies the testing directory at
    train_path into the new this directory at new_dir/user_input, then runs the
    running_stan.predict_out_of_sample function to write samples in
    new_dir/oos_samples.

    The trained output is stored in the new_dir/trained_samples folder along
    with the user input required to generate the trained samples.

    """
    click.echo(do_predict(data_path))


def do_predict(data_path: str):
    """Generate MCMC samples given a Maud output folder at train_path.

    This function creates a new directory in output_dir with a name starting
    with "maud-predict-output". It first copies the testing directory at
    train_path into the new this directory at new_dir/user_input, then runs the
    running_stan.predict_out_of_sample function to write samples in
    new_dir/oos_samples.

    The trained output is stored in the new_dir/trained_samples folder along
    with the user input required to generate the trained samples.

    """
    idata_train = az.from_json(os.path.join(data_path, "idata.json"))
    mi = load_maud_input(os.path.join(data_path, "user_input"))
    now = datetime.now().strftime("%Y%m%d%H%M%S")
    output_name = f"maud-predict_output-{mi.config.name}-{now}"
    output_path = os.path.join(data_path, output_name)
    test_samples_path = os.path.join(output_path, "test_samples")
    print("Creating output directory: " + output_path)
    os.mkdir(output_path)
    os.mkdir(test_samples_path)
    idata_predict = predict(mi, output_path, idata_train)
    idata_predict.to_json(os.path.join(output_path, "idata_predict.json"))


@cli.command("simulate")
@click.option("--output_dir", default=".", help="Where to save the output")
@click.option("-n", default=1, type=int, help="Number of simulations")
@click.argument(
    "data_path",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
)
def simulate_command(data_path, output_dir, n):
    """Generate draws from the initial values."""
    click.echo(do_simulate(data_path, output_dir, n))


@cli.command("optimize")
@click.option("--output_dir", default=".", help="Where to save the output")
@click.argument(
    "data_path",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
)
def optimize_command(data_path, output_dir):
    """Optimize the model parameters."""
    click.echo(do_optimize(data_path, output_dir))


def do_simulate(data_path, output_dir, n):
    """Generate draws from the initial values."""
    mi = load_maud_input(data_path=data_path)
    now = datetime.now().strftime("%Y%m%d%H%M%S")
    output_name = f"maud_output_sim-{mi.config.name}-{now}"
    output_path = os.path.join(output_dir, output_name)
    samples_path = os.path.join(output_path, "samples")
    ui_dir = os.path.join(output_path, "user_input")
    print("Creating output directory: " + output_path)
    os.mkdir(output_path)
    os.mkdir(samples_path)
    print(f"Copying user input from {data_path} to {ui_dir}")
    shutil.copytree(data_path, ui_dir)
    stanfit = simulate(mi, samples_path, n)
    idata = get_idata(stanfit.runset.csv_files, mi, "train")
    idata.to_json(os.path.join(output_path, "idata.json"))
    print("\n\nSimulated concentrations:")
    print(
        idata.posterior["conc_train"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
        .T
    )
    print("\n\nSimulated fluxes:")
    print(
        idata.posterior["flux_train"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
        .T
    )
    print("\n\nSimulated enzyme concentrations:")
    print(
        idata.posterior["conc_enzyme_train"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
        .T
    )
    print("\n\nSimulated reaction delta Gs:")
    print(idata.posterior["dgr_train"].mean(dim=["chain", "draw"]).to_series())
    print("\n\nSimulated concentration control coefficients")
    print(
        idata.posterior["concentration_control_matrix"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
    )
    print("\n\nSimulated flux control coefficients")
    print(
        idata.posterior["flux_control_matrix"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
    )
    print("\n\nSimulated response coefficients")
    print(
        idata.posterior["flux_response_coefficient"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
    )
    print("\n\nSimulated response coefficients")
    print(
        idata.posterior["concentration_response_coefficient"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
    )
    print("\n\nSimulated measurements:")
    for var in ["yrep_conc_train", "yrep_flux_train"]:
        if var in idata.posterior_predictive.data_vars:
            print(
                idata.posterior_predictive[var]
                .mean(dim=["chain", "draw"])
                .to_series()
            )
    print("\n\nSimulated log likelihoods:")
    for var in ["llik_conc_train", "llik_flux_train"]:
        if var in idata.log_likelihood.data_vars:
            print(
                idata.log_likelihood[var]
                .mean(dim=["chain", "draw"])
                .to_series()
            )
    print("\n\nSimulated allostery terms:")
    print(
        idata.posterior["allostery_train"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
        .T
    )
    print("\n\nSimulated reversibility terms:")
    print(
        idata.posterior["reversibility_train"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
        .T
    )
    print("\n\nSimulated saturation terms:")
    print(
        idata.posterior["saturation_train"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
        .T
    )
    if mi.kinetic_model.phosphorylations is not None:
        print("\n\nSimulated phosphorylation terms:")
        print(
            idata.posterior["phosphorylation_train"]
            .mean(dim=["chain", "draw"])
            .to_series()
            .unstack()
            .T
        )
    print("\n\nSimulated membrane potential:")
    print(
        idata.posterior["psi_train"].mean(dim=["chain", "draw"]).to_series().T
    )
    print("\n\nSimulated steady state deviation:")
    # TODO: put the dimensions directly in the idata object
    steady_dev = (
        idata.posterior["steady_dev"].mean(dim=["chain", "draw"]).to_series().T
    )
    steady_dev.index = steady_dev.index.set_names("experiment", level=0)
    steady_dev.index = steady_dev.index.set_names("metabolite", level=1)
    balanced_metabolites = [m.id for m in mi.kinetic_model.mics if m.balanced]
    steady_dev.index = steady_dev.index.set_levels(
        [
            idata.posterior.experiments,
            [
                m
                for m in idata.posterior.mics.to_series()
                if str(m) in balanced_metabolites
            ],
        ]
    )
    print(steady_dev)
    return output_path


def do_optimize(data_path, output_dir):
    """Optimize the model parameters."""
    mi = load_maud_input(data_path=data_path)
    now = datetime.now().strftime("%Y%m%d%H%M%S")
    output_name = f"maud_output_opt-{mi.config.name}-{now}"
    output_path = os.path.join(output_dir, output_name)
    samples_path = os.path.join(output_path, "samples")
    ui_dir = os.path.join(output_path, "user_input")
    print("Creating output directory: " + output_path)
    os.mkdir(output_path)
    os.mkdir(samples_path)
    print(f"Copying user input from {data_path} to {ui_dir}")
    shutil.copytree(data_path, ui_dir)
    stanfit = optimize(mi, samples_path)
    idata = get_idata(stanfit.runset.csv_files, mi, "train")
    idata.to_json(os.path.join(output_path, "idata.json"))
    print("\n\nSimulated concentrations:")
    print(
        idata.posterior["conc_train"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
        .T
    )
    print("\n\nSimulated fluxes:")
    print(
        idata.posterior["flux_train"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
        .T
    )
    print("\n\nSimulated kms:")
    print(idata.posterior["km"].mean(dim=["chain", "draw"]).to_series())
    print("\n\nSimulated enzyme concentrations:")
    print(
        idata.posterior["conc_enzyme_train"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
        .T
    )
    print("\n\nSimulated reaction delta Gs:")
    print(idata.posterior["dgr_train"].mean(dim=["chain", "draw"]).to_series())
    print("\n\nSimulated measurements:")
    print(
        idata.posterior_predictive["yrep_conc_train"]
        .mean(dim=["chain", "draw"])
        .to_series()
    )
    print("\n\nSimulated measurements:")
    for var in ["yrep_conc_train", "yrep_flux_train"]:
        if var in idata.posterior_predictive.data_vars:
            print(
                idata.posterior_predictive[var]
                .mean(dim=["chain", "draw"])
                .to_series()
            )
    print("\n\nSimulated log likelihoods:")
    for var in ["llik_conc_train", "llik_flux_train"]:
        if var in idata.log_likelihood.data_vars:
            print(
                idata.log_likelihood[var]
                .mean(dim=["chain", "draw"])
                .to_series()
            )
    print("\n\nSimulated allostery terms:")
    print(
        idata.posterior["allostery_train"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
        .T
    )
    print("\n\nSimulated reversibility terms:")
    print(
        idata.posterior["reversibility_train"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
        .T
    )
    print("\n\nSimulated saturation terms:")
    print(
        idata.posterior["saturation_train"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
        .T
    )
    if mi.kinetic_model.phosphorylations is not None:
        print("\n\nSimulated phosphorylation terms:")
        print(
            idata.posterior["phosphorylation_train"]
            .mean(dim=["chain", "draw"])
            .to_series()
            .unstack()
            .T
        )
    print("\n\nSimulated membrane potential:")
    print(
        idata.posterior["psi_train"].mean(dim=["chain", "draw"]).to_series().T
    )
    print("\n\nSimulated steady state deviation:")
    # TODO: put the dimensions directly in the idata object
    steady_dev = (
        idata.posterior["steady_dev"].mean(dim=["chain", "draw"]).to_series().T
    )
    steady_dev.index = steady_dev.index.set_names("experiment", level=0)
    steady_dev.index = steady_dev.index.set_names("metabolite", level=1)
    balanced_metabolites = [m.id for m in mi.kinetic_model.mics if m.balanced]
    steady_dev.index = steady_dev.index.set_levels(
        [
            idata.posterior.experiments,
            [
                m
                for m in idata.posterior.mics.to_series()
                if str(m) in balanced_metabolites
            ],
        ]
    )
    print(steady_dev)
    return output_path


@cli.command("variational")
@click.option("--output_dir", default=".", help="Where to save Maud's output")
@click.argument(
    "data_path",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
)
def variational_command(data_path, output_dir):
    """Generate MCMC samples given a user input directory.

    This function creates a new directory in output_dir with a name starting
    with "maud_output". It first copies the directory at data_path into the new
    this directory at new_dir/user_input, then runs the running_stan.sample
    function to write samples in new_dir/samples. Finally it prints the results
    of cmdstanpy's diagnose and summary methods.

    """
    click.echo(do_variational(data_path, output_dir))


def do_variational(data_path, output_dir):
    """Generate MCMC samples given a user input directory.

    This function creates a new directory in output_dir with a name starting
    with "maud_output". It first copies the directory at data_path into the new
    this directory at new_dir/user_input, then runs the running_stan.sample
    function to write samples in new_dir/samples. Finally it prints the results
    of cmdstanpy's diagnose and summary methods.

    """
    mi = load_maud_input(data_path)
    now = datetime.now().strftime("%Y%m%d%H%M%S")
    output_name = f"maud_output_vi-{mi.config.name}-{now}"
    output_path = os.path.join(output_dir, output_name)
    samples_path = os.path.join(output_path, "samples")
    ui_dir = os.path.join(output_path, "user_input")
    print("Creating output directory: " + output_path)
    os.mkdir(output_path)
    os.mkdir(samples_path)
    print(f"Copying user input from {data_path} to {ui_dir}")
    shutil.copytree(data_path, ui_dir)
    variational(mi, samples_path)
    return output_path


@cli.command("pathfinder")
@click.option("--output_dir", default=".", help="Where to save Maud's output")
@click.argument(
    "data_path",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
)
def pathfinder_command(data_path, output_dir):
    """Generate draws using the pathfinder algorithm."""
    click.echo(do_pathfinder(data_path, output_dir))


def do_pathfinder(data_path, output_dir):
    """Generate draws using the pathfinder algorithm."""
    mi = load_maud_input(data_path)
    now = datetime.now().strftime("%Y%m%d%H%M%S")
    output_name = f"maud_output_pf-{mi.config.name}-{now}"
    output_path = os.path.join(output_dir, output_name)
    samples_path = os.path.join(output_path, "samples")
    ui_dir = os.path.join(output_path, "user_input")
    print("Creating output directory: " + output_path)
    os.mkdir(output_path)
    os.mkdir(samples_path)
    print(f"Copying user input from {data_path} to {ui_dir}")
    shutil.copytree(data_path, ui_dir)
    pf = pathfinder(mi, samples_path)
    inits_pf = pf.create_inits()
    for i, inits_dict in enumerate(inits_pf):
        write_stan_json(
            os.path.join(samples_path, f"inits_pathfinder-{str(i)}.json"),
            inits_dict,
        )
    return output_path


@cli.command("laplace")
@click.option("--output_dir", default=".", help="Where to save the output")
@click.argument(
    "data_path",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
)
def laplace_command(data_path, output_dir):
    """Generate approximate posterior draws using the Laplace method."""
    click.echo(do_laplace(data_path, output_dir))


def do_laplace(data_path, output_dir):
    """Generate approximate posterior draws using the Laplace method."""
    mi = load_maud_input(data_path=data_path)
    now = datetime.now().strftime("%Y%m%d%H%M%S")
    output_name = f"maud_output_laplace-{mi.config.name}-{now}"
    output_path = os.path.join(output_dir, output_name)
    samples_path = os.path.join(output_path, "samples")
    ui_dir = os.path.join(output_path, "user_input")
    print("Creating output directory: " + output_path)
    os.mkdir(output_path)
    os.mkdir(samples_path)
    print(f"Copying user input from {data_path} to {ui_dir}")
    shutil.copytree(data_path, ui_dir)
    laplace_fit = laplace(mi, samples_path)
    idata = get_idata(laplace_fit._runset.csv_files, mi, "train")
    idata.to_json(os.path.join(output_path, "idata.json"))
    print("\n\nSimulated concentrations:")
    print(
        idata.posterior["conc_train"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
        .T
    )
    print("\n\nSimulated fluxes:")
    print(
        idata.posterior["flux_train"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
        .T
    )
    print("\n\nSimulated kms:")
    print(idata.posterior["km"].mean(dim=["chain", "draw"]).to_series())
    print("\n\nSimulated enzyme concentrations:")
    print(
        idata.posterior["conc_enzyme_train"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
        .T
    )
    print("\n\nSimulated reaction delta Gs:")
    print(idata.posterior["dgr_train"].mean(dim=["chain", "draw"]).to_series())
    print("\n\nSimulated measurements:")
    print(
        idata.posterior_predictive["yrep_conc_train"]
        .mean(dim=["chain", "draw"])
        .to_series()
    )
    print("\n\nSimulated measurements:")
    for var in ["yrep_conc_train", "yrep_flux_train"]:
        if var in idata.posterior_predictive.data_vars:
            print(
                idata.posterior_predictive[var]
                .mean(dim=["chain", "draw"])
                .to_series()
            )
    print("\n\nSimulated log likelihoods:")
    for var in ["llik_conc_train", "llik_flux_train"]:
        if var in idata.log_likelihood.data_vars:
            print(
                idata.log_likelihood[var]
                .mean(dim=["chain", "draw"])
                .to_series()
            )
    print("\n\nSimulated allostery terms:")
    print(
        idata.posterior["allostery_train"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
        .T
    )
    print("\n\nSimulated reversibility terms:")
    print(
        idata.posterior["reversibility_train"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
        .T
    )
    print("\n\nSimulated saturation terms:")
    print(
        idata.posterior["saturation_train"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
        .T
    )
    if mi.kinetic_model.phosphorylations is not None:
        print("\n\nSimulated phosphorylation terms:")
        print(
            idata.posterior["phosphorylation_train"]
            .mean(dim=["chain", "draw"])
            .to_series()
            .unstack()
            .T
        )
    print("\n\nSimulated membrane potential:")
    print(
        idata.posterior["psi_train"].mean(dim=["chain", "draw"]).to_series().T
    )
    print("\n\nSimulated steady state deviation:")
    # TODO: put the dimensions directly in the idata object
    steady_dev = (
        idata.posterior["steady_dev"].mean(dim=["chain", "draw"]).to_series().T
    )
    steady_dev.index = steady_dev.index.set_names("experiment", level=0)
    steady_dev.index = steady_dev.index.set_names("metabolite", level=1)
    balanced_metabolites = [m.id for m in mi.kinetic_model.mics if m.balanced]
    steady_dev.index = steady_dev.index.set_levels(
        [
            idata.posterior.experiments,
            [
                m
                for m in idata.posterior.mics.to_series()
                if str(m) in balanced_metabolites
            ],
        ]
    )
    print(steady_dev)
    return output_path
