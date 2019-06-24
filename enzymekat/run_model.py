import arviz
import numpy as np
import os
import pandas as pd
from cmdstanpy import compile_model, sample, jsondump
from python_modules import enzymekat_data, code_generation_commands
from python_modules.conversion import sem_pct_to_lognormal_sigma

MODEL_NAME = 'toy'
REL_TOL = 1e-13
ABS_TOL = 1e-9
MAX_STEPS = int(1e9)
LIKELIHOOD = 0
N_SAMPLES = 300
N_WARMUP = 300
N_CHAINS = 4
N_CORES = 4
REFRESH = 10
STEADY_STATE_TIME = 100

RELATIVE_PATHS = {
    'data': f'../data/in/{MODEL_NAME}_data.toml',
    'stan_includes': 'stan_code',
    'stan_model': f'stan_code/inference_model_{MODEL_NAME}.stan',
    'input_data_file': f'../data/stan_records/input_data_{MODEL_NAME}.json',
    'output_data': f'../data/model_output{MODEL_NAME}.csv',
    'output_infd': f'../data/infd_{MODEL_NAME}.nc',
}


def match_string_to_file(s: str, path: str):
    return open(path, 'r').read() == s

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
    input_data = {
        'N_metabolite': len(metabolite_names),
        'N_kinetic_parameter': len(kinetic_parameters),
        'N_thermodynamic_parameter': len(thermodynamic_parameters),
        'N_reaction': len(reaction_names),
        'N_known_real': len(data.known_reals),
        'N_experiment': len(data.experiments),
        'N_measurement': len(data.measurements),
        'metabolite_ix': data.measurements['metabolite_code'].values.tolist(),
        'reaction_ix': data.measurements['reaction_code'].values.tolist(),
        'experiment_ix': data.measurements['experiment_code'].values.tolist(),
        'is_flux': data.measurements['type'].eq('flux').astype(int).values.tolist(),
        'measurement': data.measurements['value'].values.tolist(),
        'measurement_scale': data.measurements['scale'].values.tolist(),
        'known_reals': data.known_reals.drop('stan_code', axis=1).values.tolist(),
        'prior_location_kinetic': kinetic_parameters['loc'].values.tolist(),
        'prior_scale_kinetic': kinetic_parameters['scale'].values.tolist(),
        'prior_location_thermodynamic': thermodynamic_parameters['loc'].values.tolist(),
        'prior_scale_thermodynamic': thermodynamic_parameters['scale'].values.tolist(),
        'initial_concentration': initial_concentration.values.tolist(),
        'initial_time': 0,
        'steady_time': STEADY_STATE_TIME,
        'rel_tol': REL_TOL,
        'abs_tol': ABS_TOL,
        'max_steps': MAX_STEPS,
        'LIKELIHOOD': LIKELIHOOD
    }

    # dump input data
    jsondump(paths['input_data_file'], input_data)

    # compile model if necessary
    stan_code = code_generation_commands.create_stan_model(data)
    if not match_string_to_file(stan_code, paths['stan_model']):
        with open(paths['stan_model'], 'w') as f:
            f.write(stan_code)
        model = compile_model(
            paths['stan_model'],
            include_paths=[paths['stan_includes']],
            overwrite=True
        )
    else:
        model = compile_model(
            paths['stan_model'],
            include_paths=[paths['stan_includes']]
        )
        

    # run model
    posterior_samples = sample(
        model,
        chains=N_CHAINS,
        cores=4,
        data_file=paths['input_data_file'],
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
            'kinetic_parameter_names': list(kinetic_parameters['name']),
            'metabolite_names': metabolite_names
        },
        dims={
            'kinetic_parameters': ['kinetic_parameter_names'],
            'measurement_hat': ['metabolite_names']
        })
    infd.to_netcdf(paths['output_infd'])
    print(arviz.summary(infd))
