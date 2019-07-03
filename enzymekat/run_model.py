import arviz
import numpy as np
import os
import pandas as pd
from cmdstanpy import compile_model, sample, jsondump, summary
from python_modules import enzymekat_data, code_generation_commands, utils
from python_modules.conversion import sem_pct_to_lognormal_sigma

MODEL_NAME = 'linear'
REL_TOL = 1e-13
ABS_TOL = 1e-12
MAX_STEPS = int(1e9)
LIKELIHOOD = 1
N_SAMPLES = 40
N_WARMUP = 40
N_CHAINS = 4
N_CORES = 4
REFRESH = 10
STEADY_STATE_TIME = 200

RELATIVE_PATHS = {
    'data': f'../data/in/{MODEL_NAME}_data.toml',
    'stan_includes': 'stan_code',
    'stan_model': f'stan_code/autogen/inference_model_{MODEL_NAME}.stan',
    'input_data_file': f'../data/stan_records/input_data_{MODEL_NAME}.json',
    'output_data': f'../data/out/model_output_{MODEL_NAME}.csv',
    'output_infd': f'../data/infd_{MODEL_NAME}.nc',
}


if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    paths = {k: os.path.join(here, v) for k, v in RELATIVE_PATHS.items()}

    # define input data
    data = enzymekat_data.from_toml(paths['data'])
    metabolite_names = data.stoichiometry.columns
    reaction_names = data.stoichiometry.index
    initial_concentration = pd.Series({m: 1 for m in metabolite_names})
    input_data = {
        'N_metabolite': len(metabolite_names),
        'N_kinetic_parameter': len(data.kinetic_parameters),
        'N_thermodynamic_parameter': len(data.thermodynamic_parameters),
        'N_reaction': len(reaction_names),
        'N_known_real': len(data.known_reals),
        'N_experiment': len(data.experiments),
        'N_measurement_conc': len(data.concentration_measurements),
        'N_measurement_flux': len(data.flux_measurements),
        'metabolite_ix': data.concentration_measurements['metabolite_code'].values.tolist(),
        'experiment_ix_conc': data.concentration_measurements['experiment_code'].values.tolist(),
        'measurement_conc': data.concentration_measurements['value'].values.tolist(),
        'measurement_scale_conc': data.concentration_measurements['scale'].values.tolist(),
        'reaction_ix': data.flux_measurements['reaction_code'].values.tolist(),
        'experiment_ix_flux': data.flux_measurements['experiment_code'].values.tolist(),
        'measurement_flux': data.flux_measurements['value'].values.tolist(),
        'measurement_scale_flux': data.flux_measurements['scale'].values.tolist(),
        'known_reals': data.known_reals.drop('stan_code', axis=1).values.tolist(),
        'prior_location_kinetic': data.kinetic_parameters['loc'].values.tolist(),
        'prior_scale_kinetic': data.kinetic_parameters['scale'].values.tolist(),
        'prior_location_thermodynamic': data.thermodynamic_parameters['loc'].values.tolist(),
        'prior_scale_thermodynamic': data.thermodynamic_parameters['scale'].values.tolist(),
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
    if not utils.match_string_to_file(stan_code, paths['stan_model']):
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
    inits = {
        'kinetic_parameters': np.exp(data.kinetic_parameters['loc']).tolist(),
        'thermodynamic_parameters': data.thermodynamic_parameters['loc'].tolist()
    }
    posterior_samples = sample(
        model,
        data=paths['input_data_file'],
        chains=N_CHAINS,
        cores=4,
        inits=inits,
        show_progress=True,
        csv_output_file=paths['output_data'],
        sampling_iters=N_SAMPLES,
        warmup_iters=N_WARMUP,
        max_treedepth=15
    )
    print(summary(posterior_samples))
