""" Simulate some measurements based on known parameter values """

import arviz
import numpy as np
import os
import pandas as pd
import cmdstanpy
from enzymekat import code_generation, data_model, utils

RELATIVE_PATHS = {
    'stan_includes': 'stan_code',
    'stan_autogen': 'stan_code/autogen',
    'stan_records': '../data/stan_records',
    'data_out': '../data/out',
}

def simulate(
        data_path: str,
        steady_state_time: float,
        rel_tol: float,
        abs_tol: float,
        max_steps: int
):
    model_name = os.path.splitext(os.path.basename(data_path))[0]
    here = os.path.dirname(os.path.abspath(__file__))
    paths = {k: os.path.join(here, v) for k, v in RELATIVE_PATHS.items()}

    # define input data
    data = data_model.from_toml(data_path)
    metabolite_names = data.stoichiometry.columns
    reaction_names = data.stoichiometry.index
    initial_concentration = pd.Series({m: 1 for m in metabolite_names})
    input_data = {
        'N_metabolite': len(metabolite_names),
        'N_param': len(data.parameters),
        'N_reaction': len(reaction_names),
        'N_experiment': len(data.experiments),
        'N_known_real': len(data.known_reals),
        'N_measurement_flux': len(data.flux_measurements),
        'N_measurement_conc': len(data.concentration_measurements),
        'metabolite_ix': data.concentration_measurements['metabolite_code'].values.tolist(),
        'experiment_ix_conc': data.concentration_measurements['experiment_code'].values.tolist(),
        'measurement_scale_conc': data.concentration_measurements['scale'].values.tolist(),
        'reaction_ix': data.flux_measurements['reaction_code'].values.tolist(),
        'experiment_ix_flux': data.flux_measurements['experiment_code'].values.tolist(),
        'measurement_scale_flux': data.flux_measurements['scale'].values.tolist(),
        'known_reals': data.known_reals.drop('stan_code', axis=1).values.tolist(),
        'params': data.parameters['true_value'].tolist(),
        'initial_concentration': initial_concentration.values.tolist(),
        'initial_time': 0,
        'steady_time': steady_state_time,
        'rel_tol': rel_tol,
        'abs_tol': abs_tol,
        'max_steps': max_steps,
    }

    # dump input data
    input_file = os.path.join(
        paths['stan_records'], f'simulation_input_data_{model_name}.json'
    )
    cmdstanpy.jsondump(input_file, input_data)

    # compile model if necessary
    stan_code = code_generation.create_stan_model(data, template='simulation')
    stan_file = os.path.join(
        paths['stan_autogen'], f'simulation_model_{model_name}.stan'
    )
    if not utils.match_string_to_file(stan_code, stan_file):
        with open(stan_file, 'w') as f:
            f.write(stan_code)
        model = cmdstanpy.compile_model(
            stan_file,
            include_paths=[paths['stan_includes']],
            overwrite=True
        )
    else:
        model = cmdstanpy.compile_model(
            stan_file,
            include_paths=[paths['stan_includes']]
        )

    # get samples
    csv_output_file = os.path.join(paths['data_out'], f'simulation_output_{model_name}.csv')
    return cmdstanpy.sample(model, data=input_file, sampling_iters=1, csv_output_file=csv_output_file)
