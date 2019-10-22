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

import click

from maud import sampling


SAMPLING_DEFAULTS = {
    "f_tol_as": 1e-6,
    "rel_tol_as": 1e-9,
    "abs_tol_as": 1e-12,
    "max_steps_as": int(1e9),
    "likelihood": 1,
    "n_samples": 5,
    "n_warmup": 5,
    "n_chains": 4,
    "n_cores": 4,
    "time_step": 0.05,
}
RELATIVE_PATH_EXAMPLE = "../../data/in/linear.toml"


def get_example_path(relative_path_example):
    """Get absolute path to file containing example input."""

    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(here, relative_path_example)


@click.group()
@click.help_option("--help", "-h")
def cli():
    """Use Maud's command line interface."""


pass


@cli.command()
@click.option(
    "--f_tol",
    default=SAMPLING_DEFAULTS["f_tol_as"],
    help="Algebra solver's functional tolerance parameter",
)
@click.option(
    "--rel_tol",
    default=SAMPLING_DEFAULTS["rel_tol_as"],
    help="Algebra solver's absolute tolerance parameter",
)
@click.option(
    "--max_steps",
    default=SAMPLING_DEFAULTS["max_steps_as"],
    help="Algebra solver's maximum steps parameter",
)
@click.option(
    "--likelihood",
    default=SAMPLING_DEFAULTS["likelihood"],
    help="Whether (1) or not (0) to run the model in likelihood mode",
)
@click.option(
    "--n_samples",
    default=SAMPLING_DEFAULTS["n_samples"],
    help="Number of post-warmup posterior samples",
)
@click.option(
    "--n_warmup", default=SAMPLING_DEFAULTS["n_warmup"], help="Number of warmup samples"
)
@click.option(
    "--n_chains", default=SAMPLING_DEFAULTS["n_chains"], help="Number of MCMC chains"
)
@click.option(
    "--n_cores",
    default=SAMPLING_DEFAULTS["n_cores"],
    help="Number of chains to run in parallel",
)
@click.option(
    "--time_step",
    default=SAMPLING_DEFAULTS["time_step"],
    help="How far ahead the ode solver simulates",
)
@click.argument(
    "data_path",
    type=click.Path(exists=True, dir_okay=False),
    default=get_example_path(RELATIVE_PATH_EXAMPLE),
)
def sample(data_path, **kwargs):
    """Sample from the model defined by the data at data_path."""
    stanfit = sampling.sample(data_path, **kwargs)
    print(stanfit.summary())
