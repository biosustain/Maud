import numpy as np
import os
import stan_utils
import pandas as pd
from sbml_functions import read_sbml_file

PATH_FROM_HERE_TO_SBML_FILE = '../data_in/t_brucei.xml'
PATH_FROM_HERE_TO_PRIORS_FILE = '../data_in/t_brucei_priors.csv'
REL_TOL = 1e-13
F_TOL = 1e-6
MAX_STEPS = int(1e9)

if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    sbml_path = os.path.join(here, PATH_FROM_HERE_TO_SBML_FILE)
    sbml_model = read_sbml_file(sbml_path)
    ode_metabolites = pd.Series(sbml_model.ode_metabolites)
    known_reals = pd.Series(sbml_model.known_reals)
    kinetic_parameters = pd.Series(sbml_model.kinetic_parameters)
    priors_path = os.path.join(here, PATH_FROM_HERE_TO_PRIORS_FILE)
    priors = (
        pd.read_csv(priors_path)
        .set_index('Parameter Name')
        .reindex(kinetic_parameters.index)
        .join(kinetic_parameters.rename('copasi_value'))
    )
    priors['guess'] = priors['Mode'].where((priors['mu'].notnull()) & ~(priors.index.isin(['Vm11r'])),
                                         other=priors['copasi_value'])
    priors['sigma'] = priors['sigma'].fillna(0.1)

    data = {
        'N_ode': len(ode_metabolites),
        'N_known_real': len(known_reals),
        'N_kinetic_parameter': len(kinetic_parameters),
        'initial_guess': ode_metabolites,
        'known_reals': known_reals,
        'log_kinetic_parameters': np.log(priors['guess']),
        'rel_tol': REL_TOL,
        'f_tol': F_TOL,
        'max_steps': MAX_STEPS
    }
    stan_model = stan_utils.StanModel_cache(os.path.join(here, '../stan/test_algebra_solver.stan'))
    fit = stan_model.sampling(data=data, chains=1, iter=1, algorithm='Fixed_param')
    print(fit.get_inits())
    print(fit['measurement_hat'])
