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


def sample(
        data_path: str,
        rel_tol: float,
        abs_tol: float,
        max_steps: int,
        likelihood: int,
        n_samples: int,
        n_warmup: int,
        n_chains: int,
        n_cores: int,
        steady_state_time: float
):
    model_name = os.path.splitext(os.path.basename(data_path))[0]
    here = os.path.dirname(os.path.abspath(__file__))
    paths = {k: os.path.join(here, v) for k, v in RELATIVE_PATHS.items()}

    # define input data
    data = data_model.from_toml(data_path)
    metabolite_names = data.stoichiometry.columns
    reaction_names = data.stoichiometry.index
    initial_concentration = pd.Series({m: 1 for m in metabolite_names})
    is_log_scale = data.parameters['is_kinetic'] | data.parameters['is_allosteric']
    input_data = {
        'N_metabolite': len(metabolite_names),
        'N_param': len(data.parameters),
        'N_log_scale_param': len(data.parameters.loc[is_log_scale]),
        'N_reaction': len(reaction_names),
        'N_experiment': len(data.experiments),
        'N_known_real': len(data.known_reals),
        'N_measurement_flux': len(data.flux_measurements),
        'N_measurement_conc': len(data.concentration_measurements),
        'metabolite_ix': data.concentration_measurements['metabolite_code'].values.tolist(),
        'experiment_ix_conc': data.concentration_measurements['experiment_code'].values.tolist(),
        'measurement_conc': data.concentration_measurements['value'].values.tolist(),
        'measurement_scale_conc': data.concentration_measurements['scale'].values.tolist(),
        'reaction_ix': data.flux_measurements['reaction_code'].values.tolist(),
        'experiment_ix_flux': data.flux_measurements['experiment_code'].values.tolist(),
        'measurement_flux': data.flux_measurements['value'].values.tolist(),
        'measurement_scale_flux': data.flux_measurements['scale'].values.tolist(),
        'known_reals': data.known_reals.drop('stan_code', axis=1).values.tolist(),
        'prior_location': data.parameters['loc'].values.tolist(),
        'prior_scale': data.parameters['scale'].values.tolist(),
        'is_log_scale': is_log_scale.astype(int).values.tolist(),
        'initial_concentration': initial_concentration.values.tolist(),
        'initial_time': 0,
        'steady_time': steady_state_time,
        'rel_tol': rel_tol,
        'abs_tol': abs_tol,
        'max_steps': max_steps,
        'LIKELIHOOD': likelihood
    }

    # dump input data
    input_file = os.path.join(
        paths['stan_records'], f'input_data_{model_name}.json'
    )
    cmdstanpy.utils.jsondump(input_file, input_data)

    # compile model if necessary
    stan_code = code_generation.create_stan_model(data)
    stan_file = os.path.join(
        paths['stan_autogen'], f'inference_model_{model_name}.stan'
    )
    exe_file = stan_file[:-5]
    no_need_to_compile = (
        os.path.exists(exe_file)
        and utils.match_string_to_file(stan_code, stan_file)
    )
    if no_need_to_compile:
        model = cmdstanpy.Model(stan_file=stan_file, exe_file=exe_file)
    else:
        with open(stan_file, 'w') as f:
            f.write(stan_code)
        model = cmdstanpy.Model(stan_file)
        model.compile(include_paths=[paths['stan_includes']], overwrite=True)

    # draw samples
    csv_output_file = os.path.join(paths['data_out'], f'output_{model_name}.csv')
    inits = {'params': np.where(is_log_scale, np.exp(input_data['prior_location']), input_data['prior_location'])}
    return model.sample(
        data=input_file,
        cores=4,
        chains=n_chains,
        inits=inits,
        show_progress=True,
        csv_basename=csv_output_file,
        sampling_iters=n_samples,
        warmup_iters=n_warmup,
        max_treedepth=15,
        adapt_delta=0.9,
        save_warmup=True,
        step_size=0.01
    )
