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
import os
import shutil
from datetime import datetime

import click
import pandas as pd
import toml

from maud import sampling
from maud.analysis import load_infd
from maud.io import (
    get_all_experiment_object,
    load_maud_input,
    parse_config,
    parse_toml_kinetic_model,
)
from maud.user_templates import get_inits_from_draw, get_prior_template


RELATIVE_PATH_EXAMPLE = "../../tests/data/linear"


def get_example_path(relative_path_example):
    """Get absolute path to file containing example input."""

    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(here, relative_path_example)


@click.group()
@click.help_option("--help", "-h")
def cli():
    """Use Maud's command line interface."""
    pass


def sample(data_path, output_dir):
    """Generate MCMC samples given a user input directory.

    This function creates a new directory in output_dir with a name starting
    with "maud_output". It first copies the directory at data_path into the new
    this directory at new_dir/user_input, then runs the sampling.sample
    function to write samples in new_dir/samples. Finally it prints the results
    of cmdstanpy's diagnose and summary methods.

    """
    mi = load_maud_input(data_path, mode="sample")
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
    stanfit = sampling.sample(mi, samples_path)
    print(stanfit.diagnose())
    print(stanfit.summary())
    infd = load_infd(stanfit.runset.csv_files, mi)
    infd.to_netcdf(os.path.join(output_path, "infd.nc"))
    return output_path


@cli.command("sample")
@click.option("--output_dir", default=".", help="Where to save Maud's output")
@click.argument(
    "data_path",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
    default=get_example_path(RELATIVE_PATH_EXAMPLE),
)
def sample_command(data_path, output_dir):
    """Run the sample function as a click command."""
    click.echo(sample(data_path, output_dir))


def generate_predictions(samples_path, oos_path, output_dir):
    """Generate MCMC samples given a Maud output folder.

    This function creates a new directory in output_dir with a name starting with
    "maud-oos_output". It first copies the testing directory at oos_path into the
    new this directory at new_dir/user_input, then runs the sampling.
    generate_predictions function to write samples in new_dir/oos_samples.
    The trained output is stored in the new_dir/trained_samples folder along with
    the user input required to generate the trained samples.
    """
    csvs = [
        os.path.join(samples_path, "samples", f)
        for f in os.listdir(os.path.join(samples_path, "samples"))
        if f.endswith(".csv")
    ]
    mi_oos = load_maud_input(data_path=oos_path, mode="predict")
    mi_train = load_maud_input(
        data_path=os.path.join(samples_path, "user_input"), mode="sample"
    )
    now = datetime.now().strftime("%Y%m%d%H%M%S")
    output_name = f"maud-oos_output-{mi_oos.config.name}-{now}"
    output_path = os.path.join(output_dir, output_name)
    trained_samples_path = os.path.join(output_path, "trained_samples")
    ui_dir = os.path.join(output_path, "user_input")
    oos_samples_path = os.path.join(output_path, "oos_samples")
    print("Creating output directory: " + output_path)
    os.mkdir(output_path)
    os.mkdir(oos_samples_path)
    print(f"Copying user input from {oos_path} to {ui_dir}")
    shutil.copytree(oos_path, ui_dir)
    print(f"Copying posterior_draws from {samples_path} to {ui_dir}")
    shutil.copytree(samples_path, trained_samples_path)
    print(f"Sampling {ui_dir} using samples from {samples_path}")
    sampling.generate_predictions(mi_oos, mi_train, csvs, oos_samples_path)
    return output_path


@cli.command("generate-predictions")
@click.option("--output_dir", default=".", help="Where to save Maud's output")
@click.option(
    "--oos_path",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
    help="Out of sample predictions from same model definition",
)
@click.argument(
    "samples_path",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
)
def generate_predictions_command(samples_path, oos_path, output_dir):
    """Run the sample function as a click command."""
    click.echo(generate_predictions(samples_path, oos_path, output_dir))


def simulate(data_path, output_dir, n):
    """Generate draws from the prior mean."""

    mi = load_maud_input(data_path=data_path, mode="sample")
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
    stanfit = sampling.simulate(mi, samples_path, n)
    infd = load_infd(stanfit.runset.csv_files, mi)
    infd.to_netcdf(os.path.join(output_path, "infd.nc"))
    print("\n\nSimulated concentrations:")
    print(infd.posterior["conc"].mean(dim=["chain", "draw"]).to_series().unstack().T)
    print("\n\nSimulated fluxes:")
    print(infd.posterior["flux"].mean(dim=["chain", "draw"]).to_series().unstack().T)
    print("\n\nSimulated enzyme concentrations:")
    print(
        infd.posterior["conc_enzyme"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
        .T
    )
    print("\n\nSimulated reaction delta Gs:")
    print(infd.posterior["dgrs"].mean(dim=["chain", "draw"]).to_series())
    print("\n\nSimulated measurements:")
    print(infd.posterior["yconc_sim"].mean(dim=["chain", "draw"]).to_series())
    print(infd.posterior["yflux_sim"].mean(dim=["chain", "draw"]).to_series())
    print("\n\nSimulated log likelihoods:")
    print(infd.posterior["log_lik_conc"].mean(dim=["chain", "draw"]).to_series())
    print(infd.posterior["log_lik_flux"].mean(dim=["chain", "draw"]).to_series())
    print("\n\nSimulated allostery:")
    print(
        infd.posterior["allostery"].mean(dim=["chain", "draw"]).to_series().unstack().T
    )
    print("\n\nSimulated reversibilities:")
    print(
        infd.posterior["reversibility"]
        .mean(dim=["chain", "draw"])
        .to_series()
        .unstack()
        .T
    )
    print("\n\nSimulated saturation:")
    print(
        infd.posterior["saturation"].mean(dim=["chain", "draw"]).to_series().unstack().T
    )
    # TODO(jorge): PMF posteriors should say the experiment they relate to
    print("\n\nSimulated proton-motive force:")
    print(infd.posterior["pmf"].mean(dim=["chain", "draw"]).to_series().T)
    return output_path


@cli.command("simulate")
@click.option("--output_dir", default=".", help="Where to save the output")
@click.option("-n", default=1, type=int, help="Number of simulations")
@click.argument(
    "data_path",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
    default=get_example_path(RELATIVE_PATH_EXAMPLE),
)
def simulate_command(data_path, output_dir, n):
    """Run the simulate function as a click command."""
    click.echo(simulate(data_path, output_dir, n))


def generate_prior_template(data_path):
    """Generate template for prior definitions.

    :params data_path: a path to a maud input folder with a kinetic model
    and optionally experimental input file.
    """
    config = parse_config(toml.load(os.path.join(data_path, "config.toml")))
    kinetic_model_path = os.path.join(data_path, config.kinetic_model_file)
    kinetic_model = parse_toml_kinetic_model(toml.load(kinetic_model_path))
    measurements_path = os.path.join(data_path, config.measurements_file)
    biological_config_path = os.path.join(data_path, config.biological_config_file)
    all_experiments = get_all_experiment_object(toml.load(biological_config_path))
    raw_measurements = pd.read_csv(measurements_path)
    output_name = "prior_template.csv"
    output_path = os.path.join(data_path, output_name)
    print("Creating template")
    prior_dataframe = get_prior_template(
        km=kinetic_model,
        raw_measurements=raw_measurements,
        experiments=all_experiments,
        mode="sample",
    )
    print(f"Saving template to: {output_path}")
    prior_dataframe.to_csv(output_path)
    return "Successfully generated prior template"


@cli.command("generate-prior-template")
@click.argument(
    "data_path",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
    default=get_example_path(RELATIVE_PATH_EXAMPLE),
)
def generate_prior_template_command(data_path):
    """Run the generate_prior_template function as a click command."""
    click.echo(generate_prior_template(data_path))


def generate_inits(data_path, chain, draw, warmup):
    """Generate template for init definitions.

    :params data_path: a path to a maud output folder with both samples
    and user_input folders
    :params chain: the sampling chain of the stan sampler you want to
    export
    :params draw: the sampling draw of the sampling chain you want to
    export from the start of the sampling or warmup phase
    :params warmup: indicator variable of if it is for the warmup
    or sampling phase
    """

    csvs = [
        os.path.join(data_path, "samples", f)
        for f in os.listdir(os.path.join(data_path, "samples"))
        if f.endswith(".csv")
    ]
    mi = load_maud_input(os.path.join(data_path, "user_input"), mode="sample")
    infd = load_infd(csvs, mi)
    output_name = "generated_inits.csv"
    output_path = os.path.join(data_path, output_name)
    print("Creating init")
    init_dataframe = get_inits_from_draw(infd, mi, chain, draw, warmup)
    print(f"Saving inits to: {output_path}")
    init_dataframe.to_csv(output_path)
    return "Successfully generated prior template"


@cli.command("generate-inits")
@click.argument(
    "data_path",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
)
@click.option("--chain", default=0, help="Sampling chain using python indexing")
@click.option(
    "--draw",
    default=0,
    help="Sampling draw using python indexing from start of phase",
)
@click.option("--warmup", default=0, help="0 if in sampling, 1 if in warmup phase")
def generate_inits_command(data_path, chain, draw, warmup):
    """Run the generate_inits function as a click command."""
    click.echo(generate_inits(data_path, chain, draw, warmup))


def variational(data_path, output_dir):
    """Generate MCMC samples given a user input directory.

    This function creates a new directory in output_dir with a name starting
    with "maud_output". It first copies the directory at data_path into the new
    this directory at new_dir/user_input, then runs the sampling.sample
    function to write samples in new_dir/samples. Finally it prints the results
    of cmdstanpy's diagnose and summary methods.

    """
    mi = load_maud_input(data_path, mode="sample")
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
    sampling.variational(mi, samples_path)
    return output_path


@cli.command("variational")
@click.option("--output_dir", default=".", help="Where to save Maud's output")
@click.argument(
    "data_path",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
    default=get_example_path(RELATIVE_PATH_EXAMPLE),
)
def variational_command(data_path, output_dir):
    """Run the sample function as a click command."""
    click.echo(variational(data_path, output_dir))
