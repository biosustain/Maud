"""
    Script for running a timecourse based on an sbml file and a stan model
    generated from it.

    Running the timecourse model requires integrating a stiff ODE system, which
    isn't possible with pystan at the moment due to a bug (see
    https://github.com/stan-dev/pystan/pull/518). As a workaround, this script
    uses cmdstan instead, which works but requires some configuration.

    To run this script you will need to follow these steps:

    1. Download (or `git clone`) a recent cmdstan 
       (https://github.com/stan-dev/cmdstan/releases)
    2. Run the following commands from the root of the cmdstan directory
       - `make stan-update`
       - `make clean-all`
       - `make build`
    3. Create a directory for models e.g. `mkdir models` 
    4. Create a directory for the specific timecourse model you want to 
       run, e.g. `mkdir models/t_brucei`
    5. Configure the upper case variables at the top of this script so that they 
       point at the things they should do.

    To build a pdf cmdstan manual, run `make manual`. You can also just run `make`
    to get a broad overview of how cmdstan works.

"""

import arviz
import numpy as np
import os
import pandas as pd
import pystan
from matplotlib import pyplot as plt
from python_modules import stan_utils, sbml_functions, enzymekat_data

MODEL_NAME = 'yeast'
RELATIVE_PATH_CMDSTAN = '../cmdstan'
RELATIVE_PATH_DATA = '../data'
RELATIVE_PATH_STAN_CODE = 'stan_code'
REL_TOL = 1e-13
ABS_TOL = 1e-9
MAX_STEPS = int(1e9)

N_TIMEPOINTS = 19
FIRST_TIMEPOINT = 0
TIMEPOINT_INTERVAL = 0.0001

if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    paths = {
        'data': os.path.join(
            here, RELATIVE_PATH_DATA, f'in/{MODEL_NAME}_data.toml'
        ),
        'cmdstan': os.path.join(
            here, RELATIVE_PATH_CMDSTAN
        ),
        'stan_model': os.path.join(
            here, RELATIVE_PATH_STAN_CODE, f'timecourse_model_{MODEL_NAME}.stan'
        ),
        'stan_input': os.path.join(
            here, RELATIVE_PATH_DATA, f'stan_records/timecourse_model_input_{MODEL_NAME}.Rdump'
        ),
        'inits': os.path.join(
            here, RELATIVE_PATH_DATA, f'stan_records/timeourse_inits_{MODEL_NAME}.Rdump'
        ),
        'output_data': os.path.join(
            here, RELATIVE_PATH_DATA, f'out/timecourse_output_{MODEL_NAME}.csv'
        ),
        'output_infd': os.path.join(
            here, RELATIVE_PATH_DATA, f'out/timecourse_infd_{MODEL_NAME}.nc'
        ),
        'fig': os.path.join(
            here, RELATIVE_PATH_DATA, f'out/timecourse_{MODEL_NAME}.png'
        ),
        'flux_fig': os.path.join(
            here, RELATIVE_PATH_DATA, f'out/timecourse_fluxes_{MODEL_NAME}.png'
        )

    }
    # define input data and write to file
    data = enzymekat_data.from_toml(paths['data'])
    timepoints = np.linspace(
        FIRST_TIMEPOINT, FIRST_TIMEPOINT + (TIMEPOINT_INTERVAL * N_TIMEPOINTS),
        num=N_TIMEPOINTS, endpoint=False
    )
    print(timepoints[1:])
    stan_input = {
        'N_ode': len(data.ode_metabolites),
        'N_known_real': len(data.known_reals),
        'N_kinetic_parameter': len(data.kinetic_parameters),
        'N_thermodynamic_parameter': len(data.thermodynamic_parameters),
        'N_timepoint': N_TIMEPOINTS,
        'N_reaction': len(data.reactions),
        'kinetic_parameters': np.exp(data.kinetic_parameters['prior_location'].values),
        'thermodynamic_parameters': data.thermodynamic_parameters['prior_location'].values,
        'known_reals': data.known_reals.values,
        'initial_state': data.ode_metabolites['initial_value'].values,
        'initial_time': timepoints[0],
        'timepoints': timepoints[1:],
        'rel_tol': REL_TOL,
        'abs_tol': ABS_TOL,
        'max_steps': MAX_STEPS
    }
    pystan.misc.stan_rdump(stan_input, paths['stan_input'])
    # compile model if necessary
    if not os.path.isfile(paths['stan_model'].replace('.stan', '.hpp')):
        path_from_cmdstan_to_program = os.path.relpath(
            paths['stan_model'].replace('.stan', ''), start=paths['cmdstan']
        )
        stan_utils.compile_stan_model_with_cmdstan(path_from_cmdstan_to_program)

    # run model
    method_config = 'sample algorithm=fixed_param num_warmup=0 num_samples=1'
    stan_utils.run_compiled_cmdstan_model(
        paths['stan_model'].replace('.stan', ''),
        paths['stan_input'],
        paths['output_data'],
        method_config=method_config
    )

    # put model in Arviz object
    infd = arviz.from_cmdstan(
        [paths['output_data']],
        coords={'metabolite_names': list(data.ode_metabolites['name']),
                'reaction_names': ['influx_fbp'] + list(data.reactions['name']) + ['outflux_pep'],
                'timepoint_labels': timepoints},
        dims={'ode_metabolite_sim': ['timepoint_labels', 'metabolite_names'],
              'flux_sim': ['timepoint_labels', 'reaction_names']}
    )

    infd.to_netcdf(paths['output_infd'])

    timecourse = infd.posterior['ode_metabolite_sim'].mean(dim=['chain', 'draw']).to_series().unstack()
    flux = infd.posterior['flux_sim'].mean(dim=['chain', 'draw']).to_series().unstack()

    plt.clf()
    f, axes = plt.subplots(2, 4, sharex=True, figsize=[8, 5])
    axes = axes.ravel()
    for ax, col in zip(axes, timecourse.columns):
        ax.plot(timecourse.index, timecourse[col])
        ax.set(title=col, xlabel='Time')
        if ax in [axes[0], axes[4]]:
            ax.set_ylabel('Concentration')
    plt.savefig(paths['fig'])

    plt.clf()
    f, axes = plt.subplots(3, 4, sharex=True, figsize=[8, 5])
    axes = axes.ravel()
    for ax, col in zip(axes, flux.columns):
        ax.plot(flux.index, flux[col])
        ax.set_title(col)
        if ax in [axes[0], axes[4], axes[8]]:
            ax.set_ylabel('Flux')
        if ax in axes[8:]:
            ax.set_xlabel('Time')
        else:
            ax.set_xlabel('')
    plt.tight_layout()
    plt.savefig(paths['flux_fig'])


    print('Timecourse:\n')
    print(timecourse)
    print('Fluxes:\n')
    print(flux)
    # print(f'Writing results to {timecourse_csv_path}...')
    # timecourse.to_csv(timecourse_csv_path)
    # print('Finished!')

