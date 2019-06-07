import arviz
import numpy as np
import os
import pandas as pd
from cmdstanpy.cmds import compile_model, sample
from python_modules import stan_utils, enzymekat_data
from python_modules.conversion import sem_pct_to_lognormal_sigma

MODEL_NAME = 'yeast'
REL_TOL = 1e-13
ABS_TOL = 1e-9
MAX_STEPS = int(1e9)
LIKELIHOOD = 1
N_SAMPLES = 300
N_WARMUP = 300
N_CHAINS = 4
N_CORES = 4
REFRESH = 10

RELATIVE_PATHS = {
    'data': f'../data/in/{MODEL_NAME}_data.toml',
    'stan_includes': 'stan_code',
    'stan_model': f'stan_code/inference_model_{MODEL_NAME}.stan',
    'output_data': f'../data/model_output{MODEL_NAME}.csv',
    'output_infd': f'../data/infd_{MODEL_NAME}.nc',
}

if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    paths = {k: os.path.join(here, v) for k, v in RELATIVE_PATHS.items()}

    # define input data
    data = enzymekat_data.from_toml(paths['data'])
    ode_metabolites = data.ode_metabolites.copy()
    ode_metabolites['ix_stan'] = range(1, len(ode_metabolites) + 1)
    ode_metabolites['measurement_scale'] = [
        sem_pct_to_lognormal_sigma(row['sem_pct'], row['measured_value'])
        for _, row in ode_metabolites.iterrows()
    ]
    ode_fluxes = data.reactions.copy()
    ode_fluxes['ix_stan'] = range(1, len(ode_fluxes) + 1)
    measured_metabolites = ode_metabolites.dropna(subset=['measured_value'])
    measured_flux = ode_fluxes.dropna(subset=['measured_value'])
    input_data = {
        'N_ode': len(ode_metabolites),
        'N_kinetic_parameter': len(data.kinetic_parameters),
        'N_known_real': len(data.known_reals),
        'N_measurement': len(measured_metabolites),
        'N_flux_measurement': len(measured_flux),
        'N_reaction': len(data.reactions),
        'N_thermodynamic_parameter': len(data.thermodynamic_parameters),
        'measurement_ix': measured_metabolites['ix_stan'].tolist(),
        'measurement': measured_metabolites['measured_value'].tolist(),
        'flux_measurement_ix': measured_flux['ix_stan'].tolist(),
        'flux_measurment': measured_flux['measured_value'].tolist(),
        'prior_location_kinetic': data.kinetic_parameters['prior_location'].tolist(),
        'prior_scale_kinetic': data.kinetic_parameters['prior_scale'].tolist(),
        'prior_location_thermodynamic': data.thermodynamic_parameters['prior_location'].tolist(),
        'prior_scale_thermodynamic': data.thermodynamic_parameters['prior_scale'].tolist(),
        'measurement_scale': measured_metabolites['measurement_scale'].tolist(),
        'flux_measurment_scale': data.experiment_info['FLUX_MEASUREMENT_SCALE'],
        'known_reals': data.known_reals.tolist(),
        'initial_state': ode_metabolites['initial_value'].tolist(),
        'initial_time': 0,
        'steady_time': data.experiment_info['STEADY_STATE_TIME'],
        'rel_tol': REL_TOL,
        'abs_tol': ABS_TOL,
        'max_steps': MAX_STEPS,
        'LIKELIHOOD': LIKELIHOOD
    }

    # compile model if necessary
    model = compile_model(
        paths['stan_model'],
        include_paths=[paths['stan_includes']]
    )

    # run model
    posterior_samples = sample(
        model,
        chains=N_CHAINS,
        cores=4,
        data=input_data,
        csv_output_file=paths['output_data'],
        refresh=REFRESH,
        post_warmup_draws_per_chain=N_SAMPLES,
        warmup_draws_per_chain=N_WARMUP,
        nuts_max_depth=15
    )

    # put model output in Arviz format and save
    infd = arviz.from_cmdstanpy(
        posterior_samples,
        coords={
            'kinetic_parameter_names': list(data.kinetic_parameters['name']),
            'metabolite_names': list(ode_metabolites['name']),
            'measured_metabolite_names': list(measured_metabolites['name'])
        },
        dims={
            'kinetic_parameters': ['kinetic_parameter_names'],
            'measurement_pred': ['measured_metabolite_names'],
            'measurement_hat': ['metabolite_names']
        })
    infd.to_netcdf(paths['output_infd'])
    print(arviz.summary(infd))
