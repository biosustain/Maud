import arviz
import numpy as np
import os
import pandas as pd
import pystan
from python_modules import stan_utils, sbml_functions

MODEL_NAME = 'yeast'
RELATIVE_PATH_CMDSTAN = '../cmdstan'
RELATIVE_PATH_DATA = '../data'
RELATIVE_PATH_STAN_CODE = 'stan_code'
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
    paths = {
        'cmdstan': os.path.join(
            here, RELATIVE_PATH_CMDSTAN
        ),
        'priors': os.path.join(
            here, RELATIVE_PATH_DATA, f'in/{MODEL_NAME}_priors.csv'
        ),
        'measurements': os.path.join(
            here, RELATIVE_PATH_DATA, f'in/{MODEL_NAME}_measurements.csv'
        ),
        'stan_model': os.path.join(
            here, RELATIVE_PATH_STAN_CODE, f'inference_model_{MODEL_NAME}.stan'
        ),
        'input_data': os.path.join(
            here, RELATIVE_PATH_DATA, f'model_input_{MODEL_NAME}.Rdump'
        ),
        'inits': os.path.join(
            here, RELATIVE_PATH_DATA, f'in/inits_{MODEL_NAME}.Rdump'
        ),
        'output_data': os.path.join(
            here, RELATIVE_PATH_DATA, f'out/model_output_{MODEL_NAME}.csv'
        ),
        'output_infd': os.path.join(
            here, RELATIVE_PATH_DATA, f'out/infd_{MODEL_NAME}.nc'
        ),
        'sbml_file': os.path.join(
            here, RELATIVE_PATH_DATA, f'in/{MODEL_NAME}.xml'
        )
    }
    # real sbml file
    sbml_model_raw = sbml_functions.read_sbml_file(paths['sbml_file'])

    # read priors
    priors = pd.read_csv(paths['priors'], sep=';').set_index('parameter')

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
        sbml_model = sbml_functions.StanReadySbmlModel(sbml_model_raw.ode_metabolites,
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
        pd.read_csv(paths['measurements'], index_col='metabolite', squeeze=True)
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
    pystan.misc.stan_rdump(data, paths['input_data'])

    # define initial parameter values and write to file
    inits = {'kinetic_parameters': np.exp(priors['mu'].values)}
    pystan.misc.stan_rdump(inits, paths['inits'])

    # compile model if necessary
    if not os.path.isfile(paths['stan_model'].replace('.stan', '.hpp')):
        path_from_cmdstan_to_program = os.path.relpath(
            paths['stan_model'].replace('.stan', ''), start=paths['cmdstan']
        )
        stan_utils.compile_stan_model_with_cmdstan(path_from_cmdstan_to_program)

    # run model
    method_config = f"""sample algorithm=hmc engine=nuts max_depth=15 \
                        num_samples={N_SAMPLES} \
                        num_warmup={N_WARMUP}"""
    stan_utils.run_compiled_cmdstan_model(
        paths['stan_model'].replace('.stan', ''),
        paths['input_data'],
        paths['output_data'],
        method_config=method_config,
        init_config=f"init={paths['inits']}",
        refresh_config="refresh=5",
        chains=N_CHAINS
    )

    # put model output in Arviz format and save
    infd = arviz.from_cmdstan(
        [paths['output_data']],
        coords={
            'kinetic_parameter_names': list(kinetic_parameters.index),
            'measurement_names': list(ode_metabolites.index)
        },
        dims={
            'kinetic_parameters': ['kinetic_parameter_names'],
            'measurement_pred': ['measurement_names'],
            'measurement_hat': ['measurement_names']
        })
    infd.to_netcdf(paths['output_infd'])
    print(arviz.summary(infd))
