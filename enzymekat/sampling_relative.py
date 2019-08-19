""" Inferes kinetic parameters and metabolite concentrations given fluxes,
priors on kinetic parametersm, and relative metabolomics data."""

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
    ed = data_model.from_toml(data_path)
    metabolites = ed.metabolites
    metabolite_names = ed.stoichiometry.columns
    unbalanced_metabolites = ed.metabolites.query('is_unbalanced')
    balanced_metabolites = ed.metabolites.query('~is_unbalanced')
    reaction_names = ed.stoichiometry.index
    unbalanced_loc, unbalanced_scale = (
        ed.unbalanced_metabolite_priors
        .set_index(['metabolite_code', 'experiment_code'])
        [col]
        .unstack()
        .values
        .tolist()
        for col in ['loc', 'scale']
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
        'flux_measurement': ed.flux_measurements['value'].values,
        'concentration_measurement': ed.concentration_measurements['value'].values,
        'flux_measurement_scale': ed.flux_measurements['scale'].values,
        'concentration_measurement_scale': ed.concentration_measurements['scale'].values,
        'prior_location_kinetic_parameter': ed.kinetic_parameters['loc'].values,
        'prior_scale_kinetic_parameter': ed.kinetic_parameters['scale'].values,
        'prior_location_unbalanced': unbalanced_loc,
        'prior_scale_unbalanced': unbalanced_scale,
        'known_reals': ed.known_reals.drop('stan_code', axis=1).values.tolist(),
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
    stan_code = code_generation.create_stan_model(ed, template='relative')
    stan_file = os.path.join(
        paths['stan_autogen'], f'relative_model_{model_name}.stan'
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
    csv_output_file = os.path.join(paths['data_out'], f'relative_output_{model_name}.csv')

    fit = model.sample(
        data=input_file,
        cores=4,
        chains=n_chains,
        show_progress=True,
        csv_basename=csv_output_file,
        sampling_iters=n_samples,
        warmup_iters=n_warmup,
        max_treedepth=15,
        adapt_delta=0.8,
        save_warmup=True
    )

    infd_posterior = (
        arviz.from_cmdstanpy(
            posterior=fit,
            posterior_predictive=['simulated_flux_measurement', 'simulated_concentration_measurement'],
            observed_data={'simulated_flux_measurement': input_data['flux_measurement'],
                           'simulated_concentration_measurement': input_data['concentration_measurement']},
            coords={
                'reactions': reaction_names,
                'metabolites': metabolite_names,
                'experiments': [x['label'] for x in ed.experiments],
                'parameter_names':  (
                [f'{x}-{y}' for x, y in zip(ed.kinetic_parameters['reaction'], ed.kinetic_parameters['parameter'])]
                ),
                'flux_measurements': (
                    [f'{x}-{y}' for x, y in zip(ed.flux_measurements['reaction_code'], ed.flux_measurements['experiment_code'])]
                ),
                'concentration_measurements': (
                    [f'{x}-{y}' for x, y in zip(ed.concentration_measurements['metabolite_code'], ed.concentration_measurements['experiment_code'])]
                )
            },
            dims={
                'concentration': ['metabolites', 'experiments'],
                'flux': ['reactions', 'experiments'],
                'kinetic_parameter': ['parameter_names'],
                'simulated_concentration_measurement': ['concentration_measurements'],
                'simulated_flux_measurement': ['flux_measurements'],
                'scaled_concentration': ['metabolites', 'experiments']

            }
        )
    )

    infd_posterior.to_netcdf(os.path.join(
        paths['data_out'], f'relative_inference_{model_name}.nc'
    ))

    return fit
