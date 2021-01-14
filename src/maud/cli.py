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
from datetime import datetime

import os

import click

from maud import sampling
from maud.io import load_maud_input_from_toml

import toml


RELATIVE_PATH_EXAMPLE = "../../tests/data/inputs/linear"


def get_example_path(relative_path_example):
    """Get absolute path to file containing example input."""

    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(here, relative_path_example)


def read_config(path):
    config_path = os.path.join(path, "config.toml")
    return toml.load(config_path)


@click.group()
@click.help_option("--help", "-h")
def cli():
    """Use Maud's command line interface."""


pass


@cli.command()
@click.option("--output_dir", default=".", help="Where to save Maud's output")
@click.argument(
    "data_path",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
    default=get_example_path(RELATIVE_PATH_EXAMPLE),
)
def sample(data_path, output_dir):
    """Sample from the model defined by the data at data_path."""
    mi = load_maud_input_from_toml(data_path)
    now = datetime.now().strftime("%Y%m%d%H%M%S")
    output_dirname = "-".join(["maud_output", mi.config.name, now])
    output_dirpath = os.path.join(output_dir, output_dirname)
    print("Creating output directory: " + output_dirpath)
    os.mkdir(output_dirpath)
    stanfit = sampling.sample(mi, output_dirpath)
    stanfit.diagnose()
    print(stanfit.summary())


# @cli.command()
# @click.option(
#     "--rel_tol",
#     default=SAMPLING_DEFAULTS["rel_tol"],
#     help="ODE solver's relative tolerance parameter",
# )
# @click.option(
#     "--abs_tol",
#     default=SAMPLING_DEFAULTS["abs_tol"],
#     help="ODE solver's absolute tolerance parameter",
# )
# @click.option(
#     "--max_num_steps",
#     default=SAMPLING_DEFAULTS["max_num_steps"],
#     help="ODE solver's maximum steps parameter",
# )
# @click.option(
#     "--timepoint",
#     default=SAMPLING_DEFAULTS["timepoint"],
#     help="How long the ODE simulates for (Units are whatever your time units are)",
# )
# @click.option("--output_dir", default=".", help="Where to save Maud's output")
# @click.argument(
#     "data_path",
#     type=click.Path(exists=True, dir_okay=False),
#     default=get_example_path(RELATIVE_PATH_EXAMPLE),
# )
# def simulate_once(data_path, **kwargs):
#     """Generate one draw using the given initial conditions for diagnostics."""
#     stanfit = sampling.simulate_once(data_path, **kwargs)
#     print(stanfit.get_drawset(params=["conc"]).T)
