import arviz
import numpy as np
import os
import pandas as pd
from cmdstanpy.cmds import compile_model, sample
from python_modules import enzymekat_data

MODEL_NAME = 'toy'
REL_TOL = 1e-13
ABS_TOL = 1e-9
MAX_STEPS = int(1e9)
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

    # define input data
    data = enzymekat_data.from_toml(paths['data'])
    metabolite_names = list(data.stoichiometry.columns)
    reaction_names = list(data.stoichiometry.index)
    kinetic_parameters = data.parameters.query("type == 'kinetic'")
    thermodynamic_parameters = data.parameters.query("type == 'thermodynamic'")
    initial_concentration = data.measurements.query("type == 'metabolite'").groupby('metabolite_code')['value'].mean()
    timepoints = np.linspace(
        FIRST_TIMEPOINT, FIRST_TIMEPOINT + (TIMEPOINT_INTERVAL * N_TIMEPOINTS),
        num=N_TIMEPOINTS, endpoint=False
    )
    stan_input = {
        'N_metabolite': len(metabolite_names),
        'N_kinetic_parameter': len(kinetic_parameters),
        'N_thermodynamic_parameter': len(data.thermodynamic_parameters),
        'N_reaction': len(data.reactions),
        'N_known_real': len(data.known_reals),
        'initial_concentration': initial_concentration.values.tolist(),
        'kinetic_parameters': kinetic_parameters['true_value'].values.tolist(),
        'thermodynamic_parameters': data.thermodynamic_parameters['true_value'].values.tolist(),
        'known_reals': data.known_reals.values.tolist(),
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
