""" Simulate some measurements based on known parameter values """

import arviz
import numpy as np
import os
import pandas as pd
import cmdstanpy
from enzymekat import code_generation, data_model, utils
import itertools

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
    ed = data_model.from_toml(data_path)
    metabolite_names = ed.stoichiometry.columns
    reaction_names = ed.stoichiometry.index
    initial_concentration = pd.Series({m: 1 for m in metabolite_names})
    input_data = {
        'N_metabolite': len(metabolite_names),
        'N_param': len(ed.parameters),
        'N_reaction': len(reaction_names),
        'N_experiment': len(ed.experiments),
        'N_known_real': len(ed.known_reals),
        'N_measurement_flux': len(ed.flux_measurements),
        'N_measurement_conc': len(ed.concentration_measurements),
        'metabolite_ix': ed.concentration_measurements['metabolite_code'].values.tolist(),
        'experiment_ix_conc': ed.concentration_measurements['experiment_code'].values.tolist(),
        'measurement_scale_conc': ed.concentration_measurements['scale'].values.tolist(),
        'reaction_ix': ed.flux_measurements['reaction_code'].values.tolist(),
        'experiment_ix_flux': ed.flux_measurements['experiment_code'].values.tolist(),
        'measurement_scale_flux': ed.flux_measurements['scale'].values.tolist(),
        'known_reals': ed.known_reals.drop('stan_code', axis=1).values.tolist(),
        'params': ed.parameters['true_value'].tolist(),
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
    cmdstanpy.utils.jsondump(input_file, input_data)

    # compile model if necessary
    stan_code = code_generation.create_stan_model(ed, template='simulation')
    stan_file = os.path.join(
        paths['stan_autogen'], f'simulation_model_{model_name}.stan'
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

    # get samples
    csv_basename = os.path.join(paths['data_out'], f'simulation_output_{model_name}.csv')
    stanfit = model.sample(data=input_file,
                        csv_basename=csv_basename,
                        adapt_engaged=False,
                        warmup_iters=0,
                        sampling_iters=1,
                        chains=1,
                        cores=1)

    # prettify simulations
    sim = stanfit.get_drawset().iloc[0]
    flux_ixs = [i for i in sim.index if 'flux_pred' in i]
    conc_ixs = [i for i in sim.index if 'conc_pred' in i]
    met_ixs = [i for i in sim.index if 'metabolite_flux' in i]
    flux_out = pd.DataFrame({
        'experiment': ed.flux_measurements['experiment_label'].values,
        'reaction': ed.flux_measurements['flux_label'].values,
        'simulated_value': sim.loc[flux_ixs].values
    })
    conc_out = pd.DataFrame({
        'experiment': ed.concentration_measurements['experiment_label'].values,
        'metabolite': ed.concentration_measurements['metabolite_label'].values,
        'simulated_value': sim.loc[conc_ixs].values
    })
    met_out_exp, met_out_met = zip(*itertools.product(
        [e['label'] for e in ed.experiments],
        list(ed.stoichiometry.columns)
    ))
    met_out = pd.DataFrame({
        'experiment': met_out_exp,
        'metabolite': met_out_met,
        'simulated_value': sim.loc[met_ixs].values

    })
    simulations = {'flux_measurements': flux_out,
                   'concentration_measurements': conc_out,
                   'metabolite_fluxes': met_out}
    return stanfit, simulations

