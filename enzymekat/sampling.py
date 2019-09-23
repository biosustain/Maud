import arviz
import numpy as np
import os
import pandas as pd
import cmdstanpy
from enzymekat import code_generation, io, utils

RELATIVE_PATHS = {
    'stan_includes': 'stan_code',
    'stan_autogen': 'stan_code/autogen',
    'stan_records': '../data/stan_records',
    'data_out': '../data/out',
}


def sample(
        data_path: str,
        f_tol: float,
        rel_tol: float,
        max_steps: int,
        likelihood: int,
        n_samples: int,
        n_warmup: int,
        n_chains: int,
        n_cores: int,
        time_step: float
):
    model_name = os.path.splitext(os.path.basename(data_path))[0]
    here = os.path.dirname(os.path.abspath(__file__))
    paths = {k: os.path.join(here, v) for k, v in RELATIVE_PATHS.items()}

    # define input data
    eki = io.load_enzymekat_input_from_toml(data_path)
    prior_df = pd.DataFrame.from_records(
        [[p.id, p.experiment_id, p.target_id, p.location, p.scale, p.target_type]
         for p in eki.priors.values()],
        columns=['id', 'experiment_id', 'target_id', 'location', 'scale', 'target_type']
    )
    metabolites = eki.kinetic_model.metabolites
    reactions = eki.kinetic_model.reactions
    enzymes = {k: v for r in reactions.values() for k, v in r.enzymes.items()}
    balanced_metabolites = {k: v for k, v in metabolites.items() if v.balanced}
    unbalanced_metabolites = {k: v for k, v in metabolites.items() if not v.balanced}
    unbalanced_metabolite_priors, kinetic_parameter_priors, enzyme_priors = (
        prior_df.loc[lambda df: df['target_type'] == target_type]
        for target_type in ['unbalanced_metabolite', 'kinetic_parameter', 'enzyme']
    )
    prior_loc_unb, prior_loc_enzyme, prior_scale_unb, prior_scale_enzyme = (
        df.set_index(['target_id', 'experiment_id'])[col].unstack()
        for col in ['location', 'scale']
        for df in [unbalanced_metabolite_priors, enzyme_priors]
    )
    metabolite_measurements, reaction_measurements = (
        pd.DataFrame([
            [exp.id, meas.target_id, meas.value, meas.uncertainty]
            for exp in eki.experiments.values()
            for meas in exp.measurements[measurement_type].values()
        ], columns=['experiment_id', 'target_id', 'value', 'uncertainty'])
        for measurement_type in ['metabolite', 'reaction']
    )
    # stan codes
    experiment_codes = utils.codify(eki.experiments.keys())
    reaction_codes = utils.codify(reactions.keys())
    enzyme_codes = utils.codify(enzymes.keys())
    metabolite_codes = utils.codify(metabolites.keys())
    input_data = {
        'N_balanced': len(balanced_metabolites),
        'N_unbalanced': len(unbalanced_metabolites),
        'N_kinetic_parameter': len(kinetic_parameter_priors),
        'N_reaction': len(reactions),
        'N_enzyme': len(enzymes),
        'N_experiment': len(eki.experiments),
        'N_flux_measurement': len(reaction_measurements),
        'N_conc_measurement': len(metabolite_measurements),
        'experiment_yconc': (
            metabolite_measurements['experiment_id'].map(experiment_codes).values
        ),
        'metabolite_yconc': (
            metabolite_measurements['target_id'].map(metabolite_codes).values
        ),
        'yconc': metabolite_measurements['value'].values,
        'sigma_conc': metabolite_measurements['uncertainty'].values,
        'experiment_yflux': (
            reaction_measurements['experiment_id'].map(experiment_codes).values
        ),
        'reaction_yflux': (
            reaction_measurements['target_id'].map(reaction_codes).values
        ),
        'yflux': reaction_measurements['value'].values,
        'sigma_flux': reaction_measurements['uncertainty'].values,
        'prior_loc_kinetic_parameter': kinetic_parameter_priors['location'].values,
        'prior_scale_kinetic_parameter': kinetic_parameter_priors['scale'].values,
        'prior_loc_unbalanced': prior_loc_unb.values,
        'prior_scale_unbalanced': prior_scale_unb.values,
        'prior_loc_enzyme': prior_loc_enzyme.values,
        'prior_scale_enzyme':prior_scale_enzyme.values,
        'balanced_guess': [1. for m in range(len(balanced_metabolites))],
        'rel_tol': rel_tol,
        'f_tol': f_tol,
        'max_steps': max_steps,
        'LIKELIHOOD': likelihood
    }

    # dump input data
    input_file = os.path.join(
        paths['stan_records'], f'input_data_{model_name}.json'
    )
    cmdstanpy.utils.jsondump(input_file, input_data)

    # compile model if necessary
    stan_code = code_generation.create_stan_program(eki, 'inference', time_step)
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

    fit = model.sample(
        data=input_file,
        cores=4,
        chains=n_chains,
        csv_basename=csv_output_file,
        sampling_iters=n_samples,
        warmup_iters=n_warmup,
        max_treedepth=15,
        adapt_delta=0.8,
        save_warmup=True,
        inits={'kinetic_parameter': np.exp(kinetic_parameter_priors['location']).T.values,
               'unbalanced': np.exp(prior_loc_unb).T.values,
               'enzyme_concentration': np.exp(prior_loc_enzyme).T.values}
    )

    infd_posterior = (
        arviz.from_cmdstanpy(
            posterior=fit,
            posterior_predictive=['yflux_sim', 'yconc_sim'],
            observed_data={'yflux_sim': input_data['yflux'],
                           'yconc_sim': input_data['yconc']},
            coords={
                'reactions': list(reaction_codes.keys()),
                'metabolites': list(metabolite_codes.keys()),
                'experiments': list(experiment_codes.keys()),
                'enzymes': list(enzyme_codes.keys()),
                'kinetic_parameter_names': kinetic_parameter_priors['id'].tolist(),
                'reaction_measurements': [
                    str(experiment_codes[row['experiment_id']]) + '_' + str(reaction_codes[row['target_id']])
                    for _, row in reaction_measurements.iterrows()
                ],
                'metabolite_measurements': [
                    str(experiment_codes[row['experiment_id']]) + '_' + str(metabolite_codes[row['target_id']])
                    for _, row in metabolite_measurements.iterrows()
                ],
            },
            dims={
                'conc': ['experiments', 'metabolites'],
                'flux': ['experiments', 'reactions'],
                'kinetic_parameter': ['kinetic_parameter_names'],
                'enzyme_concentration': ['experiments', 'enzymes'],
                'yconc_sim': ['metabolite_measurements'],
                'yflux_sim': ['reaction_measurements'],

            }
        )
    )

    infd_posterior.to_netcdf(os.path.join(
        paths['data_out'], f'model_inference_{model_name}.nc'
    ))

    return fit
