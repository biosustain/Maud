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
    metabolites = ed.metabolites
    metabolite_names = ed.stoichiometry.columns
    unbalanced_metabolites = ed.metabolites.query('is_unbalanced')
    balanced_metabolites = ed.metabolites.query('~is_unbalanced')
    reaction_names = ed.stoichiometry.index
    concentration_unbalanced = (
        ed.unbalanced_metabolite_priors
        .set_index(['metabolite_code', 'experiment_code'])
        ['true_value']
        .unstack()
        .values
        .tolist()
    )
    input_data = {
        'N_balanced': len(balanced_metabolites),
        'N_unbalanced': len(unbalanced_metabolites),
        'N_kinetic_parameter': len(ed.kinetic_parameters),
        'N_reaction': len(reaction_names),
        'N_experiment': len(ed.experiments),
        'N_known_real': len(ed.known_reals),
        'N_flux_measurement': len(ed.flux_measurements),
        'N_concentration_measurement': len(ed.concentration_measurements),
        'pos_balanced': balanced_metabolites['stan_code'].values,
        'pos_unbalanced': unbalanced_metabolites['stan_code'].values,
        'ix_experiment_concentration_measurement': ed.concentration_measurements['experiment_code'].values,
        'ix_experiment_flux_measurement': ed.flux_measurements['experiment_code'].values,
        'ix_metabolite_concentration_measurement': ed.concentration_measurements['metabolite_code'].values,
        'ix_reaction_flux_measurement': ed.flux_measurements['reaction_code'].values,
        'flux_measurement_scale': ed.flux_measurements['scale'].values,
        'concentration_measurement_scale': ed.concentration_measurements['scale'].values,
        'kinetic_parameter': ed.kinetic_parameters['true_value'].values,
        'concentration_unbalanced': concentration_unbalanced,
        'known_reals': ed.known_reals.drop('stan_code', axis=1).values.tolist(),
        'initial_time': 0,
        'steady_time': steady_state_time,
        'rel_tol': rel_tol,
        'abs_tol': abs_tol,
        'max_steps': max_steps
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
    sim_conc, sim_flux, sim_bal_rocs = (
        stanfit.get_drawset([var]).iloc[0].values
        for var in ['simulated_concentration_measurement',
                    'simulated_flux_measurement',
                    'balanced_metabolite_rate_of_change']
    )
    out_conc = (
        ed.concentration_measurements[['metabolite_label', 'experiment_label']]
        .rename(columns={'metabolite_label': 'metabolite',
                         'experiment_label': 'experiment'})
        .assign(simulated_measurement=sim_conc)
        .round(2)
    )
    out_flux = (
        ed.flux_measurements[['flux_label', 'experiment_label']]
        .rename(columns={'flux_label': 'reaction',
                         'experiment_label': 'experiment'})
        .assign(simulated_measurement=sim_flux)
        .round(2)

    )
    if any(x > 0.01 for x in sim_bal_rocs):
        print('Uh oh! Some metabolites appear to not get to steady state:\n',
              sim_bal_rocs)
    simulations = {'flux_measurements': out_flux,
                   'concentration_measurements': out_conc}
    return stanfit, simulations

