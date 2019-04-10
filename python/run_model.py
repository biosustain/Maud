import numpy as np
import os
import pandas as pd
from stan_utils import StanModel_cache_code
from sbml_functions import read_sbml_file

PATH_FROM_HERE_TO_SBML_FILE = '../data_in/t_brucei.xml'
PATH_FROM_HERE_TO_PRIORS_FILE = '../data_in/t_brucei_priors.csv'
PATH_FROM_HERE_TO_MEASUREMENT_FILE = '../data_in/t_brucei_measurements.csv'
PATH_FROM_HERE_TO_MODEL_FILE = '../stan/model.stan'
PATH_FROM_ROOT_TO_STEADY_STATE_FILE = 'stan/autogen/t_brucei.stan'
INCLUDE_FILE = 'autogen'
REL_TOL = 1e-13
F_TOL = 1e-6
MAX_STEPS = int(1e9)
LIKELIHOOD = 1
SIGMA_MEASUREMENT = 0.1


if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    sbml_file = os.path.join(here, PATH_FROM_HERE_TO_SBML_FILE)
    model_path = os.path.join(here, PATH_FROM_HERE_TO_MODEL_FILE)
    priors_path = os.path.join(here, PATH_FROM_HERE_TO_PRIORS_FILE)
    measurement_path = os.path.join(here, PATH_FROM_HERE_TO_MEASUREMENT_FILE)
    sbml_model = read_sbml_file(sbml_file)
    with open(model_path, 'r') as f:
        model_code = f.read().replace('REPLACE_THIS_WORD', PATH_FROM_ROOT_TO_STEADY_STATE_FILE) 
    stan_model = StanModel_cache_code(model_code)
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
        'measurement_ix': measurements['ix_stan'],
        'measurement': measurements['measured_value'],
        'prior_location': priors['mu'],
        'prior_scale': priors['sigma'],
        'sigma_measurement': SIGMA_MEASUREMENT,
        'known_reals': known_reals,
        'initial_guess': ode_metabolites,
        'rel_tol': REL_TOL,
        'f_tol': F_TOL,
        'max_steps': MAX_STEPS,
        'LIKELIHOOD': LIKELIHOOD
    }
    fit = stan_model.sampling(
        data=data,
        chains=1,
        init=[{'kinetic_parameters': kinetic_parameters}],
        verbose=True,
        seed=5,
        iter=1
    )
