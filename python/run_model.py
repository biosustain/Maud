import arviz
import numpy as np
import os
import pandas as pd
import pystan
import stan_utils
from sbml_functions import read_sbml_file, StanReadySbmlModel
from convert_sbml_to_stan import get_stan_program

MODEL_NAME = 'yeast'
PATH_FROM_HERE_TO_CMDSTAN_HOME = f'../cmdstan'
REL_TOL = 1e-13
ABS_TOL = 1e-6
MAX_STEPS = int(1e9)
LIKELIHOOD = 1
MEASUREMENT_SCALE = 0.05
N_SAMPLES = 1000
N_WARMUP = 1000
N_CHAINS = 4

if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))

    path_from_here_to_sbml_file = f'../data_in/{MODEL_NAME}.xml'
    path_from_here_to_priors_file = f'../data_in/{MODEL_NAME}_priors.csv'
    path_from_here_to_measurement_file = f'../data_in/{MODEL_NAME}_measurements.csv'
    path_from_here_to_stan_model = f'../stan/inference_model_{MODEL_NAME}.stan'
    path_from_here_to_input_data = f'../data_in/model_input_{MODEL_NAME}.Rdump'
    path_from_here_to_inits = f'../data_in/inits_{MODEL_NAME}.Rdump'
    path_from_here_to_output_data = f'../data_out/model_output_{MODEL_NAME}.csv'
    path_from_here_to_output_infd = f'../data_out/infd_{MODEL_NAME}.nc'

    sbml_file = os.path.join(here, path_from_here_to_sbml_file)
    stan_model_path = os.path.join(here, path_from_here_to_stan_model)
    priors_path = os.path.join(here, path_from_here_to_priors_file)
    measurement_path = os.path.join(here, path_from_here_to_measurement_file)
    cmdstan_home = os.path.join(here, PATH_FROM_HERE_TO_CMDSTAN_HOME)
    input_data_path = os.path.join(here, path_from_here_to_input_data)
    init_path = os.path.join(here, path_from_here_to_inits)
    output_data_path = os.path.join(here, path_from_here_to_output_data)
    output_infd_path = os.path.join(here, path_from_here_to_output_infd)

    # real sbml file
    sbml_model_raw = read_sbml_file(sbml_file)

    # read priors
    priors = pd.read_csv(priors_path).set_index('parameter')

    # check in case some parameters don't have priors
    if set(sbml_model_raw.kinetic_parameters.keys()) != set(priors.keys()):

        no_priors = {k: v
                     for k, v in sbml_model_raw.kinetic_parameters.items()
                     if k not in priors.index}
        print(f"Treating these parameters with no priors as constants:\n{no_priors}")
        new_known_reals = {**sbml_model_raw.known_reals, **no_priors}
        new_kinetic_parameters = {k: v
                                  for k, v in sbml_model_raw.kinetic_parameters.items()
                                  if k in priors.index}
        sbml_model = StanReadySbmlModel(sbml_model_raw.ode_metabolites,
                                        new_known_reals,
                                        new_kinetic_parameters,
                                        sbml_model_raw.assignment_expressions,
                                        sbml_model_raw.function_definitions,
                                        sbml_model_raw.kinetic_expressions,
                                        sbml_model_raw.ode_expressions)
    else:
        sbml_model = sbml_model_raw.copy()

    # get model input
    ode_metabolites = pd.Series(sbml_model.ode_metabolites)
    known_reals = pd.Series(sbml_model.known_reals)
    kinetic_parameters = pd.Series(sbml_model.kinetic_parameters)
    measurements = (
        pd.read_csv(measurement_path, index_col='metabolite', squeeze=True)
        .reindex(sbml_model.ode_metabolites.keys())
        .reset_index()
        .assign(ix_stan=lambda df: range(1, len(df) + 1))
        .dropna(how='any')
    )

    # define input data and write to file
    data = {
        'N_ode': len(ode_metabolites),
        'N_measurement': len(measurements),
        'N_known_real': len(known_reals),
        'N_kinetic_parameter': len(kinetic_parameters),
        'measurement_ix': measurements['ix_stan'].values,
        'measurement': measurements['measured_value'].values,
        'prior_location': priors['mu'].values,
        'prior_scale': priors['sigma'].values,
        'measurement_scale': MEASUREMENT_SCALE,
        'known_reals': known_reals.values,
        'initial_state': ode_metabolites.values,
        'initial_time': 0,
        'steady_time': 10,
        'rel_tol': REL_TOL,
        'abs_tol': ABS_TOL,
        'max_steps': MAX_STEPS,
        'LIKELIHOOD': LIKELIHOOD
    }
    pystan.misc.stan_rdump(data, input_data_path)

    # define initial parameter values and write to file
    inits = {'kinetic_parameters': np.exp(priors['mu'].values)}
    pystan.misc.stan_rdump(inits, init_path)

    # compile model if necessary
    if not os.path.isfile(stan_model_path.replace('.stan', '.hpp')):
        path_from_cmdstan_home_to_program = os.path.relpath(
            stan_model_path.replace('.stan', ''), start=cmdstan_home
        )
        stan_utils.compile_stan_model_with_cmdstan(
            path_from_cmdstan_home_to_program
        )

    # run model
    method_config = f"""sample algorithm=hmc engine=nuts max_depth=15 \
                        num_samples={N_SAMPLES} \
                        num_warmup={N_WARMUP}"""
    stan_utils.run_compiled_cmdstan_model(
        stan_model_path.replace('.stan', ''),
        input_data_path,
        output_data_path,
        method_config=method_config,
        init_config=f"init={init_path}",
        refresh_config="refresh=5",
        chains=N_CHAINS
    )

    # put model output in Arviz format and save
    infd = arviz.from_cmdstan(
        [output_data_path],
        coords={
            'kinetic_parameter_names': list(kinetic_parameters.index),
            'measurement_names': list(ode_metabolites.index)
        },
        dims={
            'kinetic_parameters': ['kinetic_parameter_names'],
            'measurement_pred': ['measurement_names'],
            'measurement_hat': ['measurement_names']
        })
    infd.to_netcdf(output_infd_path)
    print(arviz.summary(infd))
