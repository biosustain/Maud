"""Functions that are exposed to the command line interface live here."""

import click
from cmdstanpy import summary
from enzymekat import sampling, simulation
import os

SAMPLING_DEFAULTS = {
    'rel_tol': 1e-13,
    'abs_tol': 1e-12,
    'max_steps': int(1e9),
    'likelihood': 1,
    'n_samples': 40,
    'n_warmup': 40,
    'n_chains': 4,
    'n_cores': 4,
    'steady_state_time': 200
}
RELATIVE_PATH_EXAMPLE = '../data/in/linear.toml'

def get_example_path(relative_path_example):
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(here, relative_path_example)

@click.group()
@click.help_option("--help", "-h")
def cli():
    """Main entry point for enzymeKAT's command line interface."""
pass


@cli.command()
@click.option('--rel_tol', default=SAMPLING_DEFAULTS['rel_tol'],
              help="ODE solver's relative tolerance parameter")
@click.option('--abs_tol', default=SAMPLING_DEFAULTS['abs_tol'],
              help="ODE solver's absolute tolerance parameter")
@click.option('--max_steps', default=SAMPLING_DEFAULTS['max_steps'],
              help="ODE solver's maximum steps parameter")
@click.option('--likelihood', default=SAMPLING_DEFAULTS['likelihood'],
              help="Whether (1) or not (0) to run the model in likelihood mode")
@click.option('--n_samples', default=SAMPLING_DEFAULTS['n_samples'],
              help="Number of post-warmup posterior samples")
@click.option('--n_warmup', default=SAMPLING_DEFAULTS['n_warmup'],
              help="Number of warmup samples")
@click.option('--n_chains', default=SAMPLING_DEFAULTS['n_chains'],
              help="Number of MCMC chains")
@click.option('--n_cores', default=SAMPLING_DEFAULTS['n_cores'],
              help="Number of chains to run in parallel")
@click.option('--steady_state_time', default=SAMPLING_DEFAULTS['steady_state_time'],
              help="Number of time units to simulate ODEs for")
@click.argument('data_path', type=click.Path(exists=True, dir_okay=False),
                default=get_example_path(RELATIVE_PATH_EXAMPLE))
def sample(data_path, **kwargs):
    """Sample from the model defined by the data at data_path."""
    runset = sampling.sample(data_path, **kwargs)
    print(summary(runset))


@cli.command()
@click.argument('data_path', type=click.Path(exists=True, dir_okay=False),
                default=get_example_path(RELATIVE_PATH_EXAMPLE))
@click.option('--steady_state_time', default=SAMPLING_DEFAULTS['steady_state_time'],
              help="Number of time units to simulate ODEs for")
@click.option('--rel_tol', default=SAMPLING_DEFAULTS['rel_tol'],
              help="ODE solver's relative tolerance parameter")
@click.option('--abs_tol', default=SAMPLING_DEFAULTS['abs_tol'],
              help="ODE solver's absolute tolerance parameter")
@click.option('--max_steps', default=SAMPLING_DEFAULTS['max_steps'],
              help="ODE solver's maximum steps parameter")
def simulate(data_path, **kwargs):
    """Simulate measurements given parameter values from data_path."""
    runset = simulation.simulate(data_path, **kwargs)
