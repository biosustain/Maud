import arviz
import numpy as np
import os
import pandas as pd
import pystan
import stan_utils
from sbml_functions import read_sbml_file

PATH_FROM_HERE_TO_SBML_FILE = '../data_in/t_brucei.xml'
PATH_FROM_HERE_TO_PRIORS_FILE = '../data_in/t_brucei_priors.csv'
PATH_FROM_HERE_TO_MEASUREMENT_FILE = '../data_in/t_brucei_measurements.csv'
PATH_FROM_HERE_TO_STAN_MODEL_TEMPLATE = '../stan/model.stan'
PATH_FROM_HERE_TO_STAN_MODEL = '../stan/model_t_brucei.stan'
PATH_FROM_HERE_TO_CMDSTAN_HOME = '../cmdstan'
PATH_FROM_HERE_TO_INPUT_DATA = '../data_in/model_input_t_brucei.Rdump'
PATH_FROM_HERE_TO_INITS = '../data_in/inits_t_brucei.Rdump'
PATH_FROM_HERE_TO_OUTPUT_DATA = '../data_out/model_output_t_brucei.csv'
PATH_FROM_HERE_TO_OUTPUT_INFD = '../data_out/infd_t_brucei.nc'
PATH_FROM_CMDSTAN_HOME_TO_STEADY_STATE_FILE = '../stan/autogen/t_brucei.stan'
INCLUDE_FILE = 'autogen'
REL_TOL = 1e-13
ABS_TOL = 1e-6
MAX_STEPS = int(1e9)
LIKELIHOOD = 1
MEASUREMENT_SCALE = 0.05


if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    sbml_file = os.path.join(here, PATH_FROM_HERE_TO_SBML_FILE)
    stan_template_path = os.path.join(here, PATH_FROM_HERE_TO_STAN_MODEL_TEMPLATE)
    stan_model_path = os.path.join(here, PATH_FROM_HERE_TO_STAN_MODEL)
    priors_path = os.path.join(here, PATH_FROM_HERE_TO_PRIORS_FILE)
    measurement_path = os.path.join(here, PATH_FROM_HERE_TO_MEASUREMENT_FILE)
    cmdstan_home = os.path.join(here, PATH_FROM_HERE_TO_CMDSTAN_HOME)
    input_data_path = os.path.join(here, PATH_FROM_HERE_TO_INPUT_DATA)
    init_path = os.path.join(here, PATH_FROM_HERE_TO_INITS)
    output_data_path = os.path.join(here, PATH_FROM_HERE_TO_OUTPUT_DATA)
    output_infd_path = os.path.join(here, PATH_FROM_HERE_TO_OUTPUT_INFD)

    # real sbml file
    sbml_model = read_sbml_file(sbml_file)

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
    priors = (
        pd.read_csv(priors_path)
        .set_index('Parameter Name')
        .reindex(kinetic_parameters.index)
    )
    priors['mu'] = priors['mu'].fillna(0)
    priors['sigma'] = priors['sigma'].fillna(0.1)
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
    inits = {'kinetic_parameters': np.exp(priors['mu'].fillna(0)).values}
    pystan.misc.stan_rdump(data, input_data_path)
    pystan.misc.stan_rdump(inits, init_path)
    # compile model if necessary
    with open(stan_template_path, 'r') as f:
        model_code = f.read().replace(
            'REPLACE_THIS_WORD', PATH_FROM_CMDSTAN_HOME_TO_STEADY_STATE_FILE
        ) 
        f.close()
        with open(stan_model_path, 'w') as f:
            f.write(model_code)
            f.close()
    if not os.path.isfile(stan_model_path.replace('.stan', '.hpp')):
        path_from_cmdstan_home_to_program = os.path.relpath(stan_model_path.replace('.stan', ''), start=cmdstan_home)
        stan_utils.compile_stan_model_with_cmdstan(path_from_cmdstan_home_to_program)
    # run model
    extra_config = f"warmup refresh=1 random seed=326553744"
    stan_utils.run_compiled_cmdstan_model(stan_model_path.replace('.stan', ''),
                                          input_data_path,
                                          output_data_path,
                                          method_config="sample num_samples=50 num_warmup=50",
                                          refresh_config="refresh=1")
    infd = arviz.from_cmdstan([output_data_path],
                              coords={'kinetic_parameter_names': list(kinetic_parameters.index)},
                              dims={'kinetic_parameters': ['kinetic_parameter_names']})
    infd.to_netcdf(output_infd_path)
    print(arviz.summary(infd))
