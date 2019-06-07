import arviz
import numpy as np
import os
import pandas as pd
from cmdstanpy.cmds import compile_model, sample
from matplotlib import pyplot as plt
from python_modules import enzymekat_data

MODEL_NAME = 'yeast'
REL_TOL = 1e-13
ABS_TOL = 1e-9
MAX_STEPS = int(1e9)

N_TIMEPOINTS = 20
FIRST_TIMEPOINT = 0
TIMEPOINT_INTERVAL = 0.005

RELATIVE_PATHS = {
    'data': f'../data/in/{MODEL_NAME}_data.toml',
    'stan_includes': 'stan_code',
    'stan_model': f'stan_code/timecourse_model_{MODEL_NAME}.stan',
    'output_data': f'../data/out/timecourse_model_output_{MODEL_NAME}.csv',
    'output_timecourse': f'../data/out/timecourse_output_{MODEL_NAME}.csv',
    'output_flux': f'../data/out/timecourse_flux_output_{MODEL_NAME}.csv',
    'fig': f'../data/out/timecourse_{MODEL_NAME}.png',
    'flux_fig': f'../data/out/timecourse_fluxes_{MODEL_NAME}.png'
}


if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    paths = {k: os.path.join(here, v) for k, v in RELATIVE_PATHS.items()}

    # define input data and write to file
    data = enzymekat_data.from_toml(paths['data'])
    timepoints = np.linspace(
        FIRST_TIMEPOINT, FIRST_TIMEPOINT + (TIMEPOINT_INTERVAL * N_TIMEPOINTS),
        num=N_TIMEPOINTS, endpoint=False
    )
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

    # compile model if necessary
    model = compile_model(
        paths['stan_model'],
        include_paths=[paths['stan_includes']]
    )

    # run model
    posterior_samples = sample(
        model,
        chains=1,
        cores=1,
        data=stan_input,
        csv_output_file=paths['output_data'],
        post_warmup_draws_per_chain=2,
        warmup_draws_per_chain=2
    )

    # put model in Arviz object
    infd = arviz.from_cmdstanpy(
        posterior_samples,
        coords={'metabolite_names': list(data.ode_metabolites['name']),
                'reaction_names': ['influx_fbp'] + list(data.reactions['name']) + ['outflux_pep'],
                'timepoint_labels': timepoints},
        dims={'ode_metabolite_sim': ['timepoint_labels', 'metabolite_names'],
              'flux_sim': ['timepoint_labels', 'reaction_names']}
    )

    timecourse = infd.posterior['ode_metabolite_sim'].mean(dim=['chain', 'draw']).to_series().unstack()
    flux = infd.posterior['flux_sim'].mean(dim=['chain', 'draw']).to_series().unstack()

    plt.clf()
    f, axes = plt.subplots(2, 4, sharex=True, figsize=[8, 5])
    axes = axes.ravel()
    for ax, col in zip(axes, timecourse.columns):
        ax.plot(timecourse.index, timecourse[col])
        ax.set_title(col)
        if ax in [axes[0], axes[4]]:
            ax.set_ylabel('Concentration')
        if ax in axes[4:]:
            ax.set_xlabel('Time')
    plt.tight_layout()
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

    timecourse.to_csv(paths['output_timecourse'])
    timecourse.to_csv(paths['output_flux'])

    print('Timecourse:\n')
    print(timecourse)
    print('Fluxes:\n')
    print(flux)

