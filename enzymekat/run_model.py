import arviz
import numpy as np
import os
import pandas as pd
import pystan
from python_modules import stan_utils, sbml_functions, enzymekat_data

MODEL_NAME = 'yeast'
RELATIVE_PATH_CMDSTAN = '../cmdstan'
RELATIVE_PATH_DATA = '../data'
RELATIVE_PATH_STAN_CODE = 'stan_code'
REL_TOL = 1e-13
ABS_TOL = 1e-9
MAX_STEPS = int(1e9)
LIKELIHOOD = 1
MEASUREMENT_SCALE = 0.05
N_SAMPLES = 100
N_WARMUP = 100
N_CHAINS = 4

if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    paths = {
        'data': os.path.join(
            here, RELATIVE_PATH_DATA, f'in/{MODEL_NAME}_data.toml'
        ),
        'cmdstan': os.path.join(
            here, RELATIVE_PATH_CMDSTAN
        ),
        'stan_model': os.path.join(
            here, RELATIVE_PATH_STAN_CODE, f'inference_model_{MODEL_NAME}.stan'
        ),
        'stan_input': os.path.join(
            here, RELATIVE_PATH_DATA, f'stan_records/model_input_{MODEL_NAME}.Rdump'
        ),
        'inits': os.path.join(
            here, RELATIVE_PATH_DATA, f'stan_records/inits_{MODEL_NAME}.Rdump'
        ),
        'output_data': os.path.join(
            here, RELATIVE_PATH_DATA, f'out/model_output_{MODEL_NAME}.csv'
        ),
        'output_infd': os.path.join(
            here, RELATIVE_PATH_DATA, f'out/infd_{MODEL_NAME}.nc'
        )
    }
    # define input data and write to file
    data = enzymekat_data.from_toml(paths['data'])
    ode_metabolites = (
        data.ode_metabolites
        .assign(ix_stan=lambda df: range(1, len(df) + 1))
    )
    ode_fluxes = (
        data.reactions
        .assign(ix_stan=lambda df: range(1, len(df) + 1))
    )
    measured_metabolites = ode_metabolites.dropna(subset=['measured_value'])
    measured_flux = ode_fluxes.dropna(subset=['measured_value'])

    stan_input = {
        'N_ode': len(ode_metabolites),
        'N_kinetic_parameter': len(data.kinetic_parameters),
        'N_known_real': len(data.known_reals),
        'N_measurement': len(measured_metabolites),
        'N_flux_measurement': len(measured_flux),
        'N_reaction': len(data.reactions),
        'N_thermodynamic_parameter': len(data.thermodynamic_parameters),
        'measurement_ix': measured_metabolites['ix_stan'].values,
        'measurement': measured_metabolites['measured_value'].values,
        'flux_measurement_ix': measured_flux['ix_stan'].values,
        'flux_measurment': measured_flux['measured_value'].values,
        'prior_location_kinetic': data.kinetic_parameters['prior_location'].values,
        'prior_scale_kinetic': data.kinetic_parameters['prior_scale'].values,
        'prior_location_thermodynamic': data.thermodynamic_parameters['prior_location'].values,
        'prior_scale_thermodynamic': data.thermodynamic_parameters['prior_scale'].values,
        'measurement_scale': data.experiment_info['MEASUREMENT_SCALE'],
        'flux_measurment_scale': data.experiment_info['FLUX_MEASUREMENT_SCALE'],
        'known_reals': data.known_reals.values,
        'initial_state': ode_metabolites['initial_value'].values,
        'initial_time': 0,
        'steady_time': data.experiment_info['STEADY_STATE_TIME'],
        'rel_tol': REL_TOL,
        'abs_tol': ABS_TOL,
        'max_steps': MAX_STEPS,
        'LIKELIHOOD': LIKELIHOOD
    }
    pystan.misc.stan_rdump(stan_input, paths['stan_input'])

    # define initial parameter values and write to file
    inits = {'kinetic_parameters': np.exp(data.kinetic_parameters['prior_location'].values),
             'thermodynamic_parameters': data.thermodynamic_parameters['prior_location'].values}
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
        paths['stan_input'],
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
            'kinetic_parameter_names': list(data.kinetic_parameters['name']),
            'measurement_names': list(ode_metabolites['name'])
        },
        dims={
            'kinetic_parameters': ['kinetic_parameter_names'],
            'measurement_pred': ['measurement_names'],
            'measurement_hat': ['measurement_names']
        })
    infd.to_netcdf(paths['output_infd'])
    print(arviz.summary(infd))
