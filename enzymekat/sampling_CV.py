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

def reorder_indicies(indexed_list):
    current_value_indexed = indexed_list[0]
    index = 1
    output_matrix = np.zeros(len(indexed_list)).tolist()
    for ind in np.arange(0, len(indexed_list)):
        if indexed_list[ind] == current_value_indexed:
            output_matrix[ind] = index
        else:
            index += 1
            output_matrix[ind] = index
            current_value_indexed = indexed_list[ind]
    return output_matrix


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
    # partitioning data into training and validating sets
    Holdout_experiments = [1]
    Training_experiments = [2]

    flux_holdout_list = [ed.flux_measurements['experiment_code'] == x for x in Holdout_experiments][0]
    concentration_holdout_list = [ed.concentration_measurements['experiment_code'] == x for x in Holdout_experiments][0]

    flux_training_list = [ed.flux_measurements['experiment_code'] == x for x in Training_experiments][0]
    concentration_training_list = [ed.concentration_measurements['experiment_code'] == x for x in Training_experiments][0]

    flux_holdout = ed.flux_measurements[flux_holdout_list]
    concentration_holdout = ed.concentration_measurements[concentration_holdout_list]
    flux_training = ed.flux_measurements[flux_training_list]
    concentration_training = ed.concentration_measurements[concentration_training_list]

    inter_training_flux_index = reorder_indicies(flux_training['experiment_code'].values.tolist())
    inter_training_concentration_index = reorder_indicies(concentration_training['experiment_code'].values.tolist())
    inter_holdout_flux_index = reorder_indicies(flux_holdout['experiment_code'].values.tolist())
    inter_holdout_concentration_index = reorder_indicies(concentration_holdout['experiment_code'].values.tolist())


    input_data = {
        'N_balanced': len(balanced_metabolites),
        'N_unbalanced': len(unbalanced_metabolites),
        'N_training_experiments': len(Training_experiments),
        'N_holdout_experiments': len(Holdout_experiments),
        'N_kinetic_parameter': len(ed.kinetic_parameters),
        'N_reaction': len(reaction_names),
        'N_experiment': len(ed.experiments),
        'N_known_real': len(ed.known_reals),
        'N_flux_measurement_training': len(flux_training),
        'N_flux_measurement_holdout': len(flux_holdout),
        'N_concentration_measurement_training': len(concentration_training),
        'N_concentration_measurement_holdout': len(concentration_holdout),
        'N_training_experiment': len(Training_experiments),
        'N_holdout_experiments': len(Holdout_experiments),
        'Training_experiments': Training_experiments,
        'Holdout_experiments': Holdout_experiments,
        'pos_balanced': balanced_metabolites['stan_code'].values,
        'pos_unbalanced': unbalanced_metabolites['stan_code'].values,
        # 'ix_experiment_concentration_measurement': ed.concentration_measurements['experiment_code'].values,
        'ix_training_flux_experiment':
        flux_training['experiment_code'].values.tolist(),
        'ix_training_concentration_experiment': concentration_training['experiment_code'].values.tolist(),
        'ix_holdout_flux_experiment': flux_holdout['experiment_code'].values.tolist(),
        'ix_holdout_concentration_experiment': concentration_holdout['experiment_code'].values.tolist(),
        'ig_training_concentration_experiment': inter_training_concentration_index,
        'ig_training_flux_experiment': inter_training_flux_index,
        'ig_holdout_concentration_experiment': inter_holdout_concentration_index,
        'ig_holdout_flux_experiment': inter_holdout_flux_index,
        'ix_training_flux_measurement':
        flux_training['reaction_code'].values.tolist(),
        'ix_training_concentration_measurement': concentration_training['metabolite_code'].values.tolist(),
        'ix_holdout_flux_measurement': flux_holdout['reaction_code'].values.tolist(),
        'ix_holdout_concentration_measurement': concentration_holdout['metabolite_code'].values.tolist(),
        'flux_measurement_training': flux_training['value'].values.tolist(),
        'concentration_measurement_training': concentration_training['value'].values.tolist(),
        'flux_measurement_holdout': flux_holdout['value'].values.tolist(),
        'concentration_measurement_holdout': concentration_holdout['value'].values.tolist(),
        'flux_measurement_scale_training': flux_training['scale'].values.tolist(),
        'concentration_measurement_scale_training': concentration_training['scale'].values.tolist(),
        'flux_measurement_scale_holdout': flux_holdout['scale'].values.tolist(),
        'concentration_measurement_scale_holdout': concentration_holdout['scale'].values.tolist(),
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
    stan_code = code_generation.create_stan_model(ed, template='validation')
    stan_file = os.path.join(
        paths['stan_autogen'], f'inference_model_validation_{model_name}.stan'
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
    csv_output_file = os.path.join(paths['data_out'], f'output_{model_name}_validation.csv')

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

            }
        )
    )

    infd_posterior.to_netcdf(os.path.join(
        paths['data_out'], f'model_inference_{model_name}.nc'
    ))

    return fit
