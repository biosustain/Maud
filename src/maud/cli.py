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
from maud.io import load_maud_input_from_toml, parse_config, parse_toml_kinetic_model
from maud.utils import get_input_template


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
    mi = load_maud_input_from_toml(data_path)
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


def simulate(data_path, output_dir, n):
    """Generate draws from the prior mean."""

    mi = load_maud_input_from_toml(data_path)
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
    print("\nSimulated concentrations and fluxes:")
    print(infd.posterior["conc"].mean(dim=["chain", "draw"]).to_series())
    print(infd.posterior["flux"].mean(dim=["chain", "draw"]).to_series())
    print(infd.posterior["conc_enzyme"].mean(dim=["chain", "draw"]).to_series())
    print("\nSimulated measurements:")
    print(infd.posterior["yconc_sim"].mean(dim=["chain", "draw"]).to_series())
    print(infd.posterior["yflux_sim"].mean(dim=["chain", "draw"]).to_series())
    print("\nSimulated log likelihoods:")
    print(infd.posterior["log_lik_conc"].mean(dim=["chain", "draw"]).to_series())
    print(infd.posterior["log_lik_flux"].mean(dim=["chain", "draw"]).to_series())
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
    """Generate draws from the prior mean."""

    config = parse_config(toml.load(os.path.join(data_path, "config.toml")))
    kinetic_model_path = os.path.join(data_path, config.kinetic_model_file)
    kinetic_model = parse_toml_kinetic_model(toml.load(kinetic_model_path))
    experiments_path = os.path.join(data_path, config.experiments_file)
    raw_measurements = pd.read_csv(experiments_path)
    output_name = f"prior_template.csv"
    output_path = os.path.join(data_path, output_name)
    print("Creating template: " + output_path)
    prior_dataframe = get_input_template(kinetic_model, raw_measurements)
    print("Saving template")
    prior_dataframe.to_csv(output_path)
    print("Generated CSV")
    return output_path


@cli.command("generate-prior-template")
@click.argument(
    "data_path",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
    default=get_example_path(RELATIVE_PATH_EXAMPLE),
)
def generate_prior_template_command(data_path):
    """Run the simulate function as a click command."""
    click.echo(generate_prior_template(data_path))
