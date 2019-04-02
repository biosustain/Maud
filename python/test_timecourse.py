import arviz
import argparse
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import pystan
import os
from stan_utils import StanModel_cache
from sbml_functions import read_sbml_file

DEFAULT_SBML_FILE = '../data_in/t_brucei.xml'
DEFAULT_KINETICS_FILE = '../stan/autogen/t_brucei.stan'
PATH_FROM_HERE_TO_TEMPLATE = '../stan/timecourse_model_template.stan'

TIME_POINTS = np.linspace(0, 20, 101)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sbml_file',
                        type=str,
                        default=DEFAULT_SBML_FILE,
                        help='Sbml file that was used to generate the kinetics')
    parser.add_argument('--kinetics_file',
                        type=str,
                        default=DEFAULT_KINETICS_FILE,
                        help='Stan file with kinetic equations')
    args = parser.parse_args()
    path_from_here_to_kinetics_file = args.kinetics_file
    path_from_here_to_sbml_file = args.sbml_file

    # compile Stan model
    here = os.path.dirname(__file__)
    template_path = os.path.join(here, PATH_FROM_HERE_TO_TEMPLATE)
    kinetics_file_path = os.path.join(here, path_from_here_to_kinetics_file)
    with open(template_path, 'r') as f:
        model_code = f.read().replace('REPLACE_THIS_WORD', kinetics_file_path)
    model =StanModel_cache(model_code)

    # parse sbml file
    m = read_sbml_file(here + path_from_here_to_sbml_file)


    data = {
        'N_ode': len(m.ode_metabolites.values()),
        'N_derived': len(m.assignment_rules.values()),
        'N_known_real': len(m.known_reals.values()),
        'P': len(m.kinetic_parameters.values()),
        'T': len(TIME_POINTS) - 1,
        'initial_metabolite_ode': list(m.ode_metabolites.values()),
        'kinetic_parameters': list(m.kinetic_parameters.values()),
        'known_reals': list(m.known_reals.values()),
        'ts': TIME_POINTS[1:],
        't0': TIME_POINTS[0]

    }
    fit = model.sampling(data=data, algorithm='Fixed_param', iter=1, chains=1)

    infd = arviz.from_pystan(posterior=fit,
                             coords={'sim_time': TIME_POINTS,
                                     'metabolite_ode': list(m.ode_metabolites.keys()),
                                     'metabolite_derived': list(m.assignment_rules.keys())},
                             dims={'ode_metabolite_sim': ['sim_time', 'metabolite_ode'],
                                   'derived_quantity_sim': ['sim_time', 'metabolite_derived']})
    ode_out = infd.posterior['ode_metabolite_sim'].mean(dim=['chain', 'draw']).to_series().unstack()
    derived_out = infd.posterior['derived_quantity_sim'].mean(dim=['chain', 'draw']).to_series().unstack()
    out = ode_out.join(derived_out)
    out = out.sort_index()
    f, axes = plt.subplots(5, 6, sharex=True, figsize=[15, 15])
    axes = axes.ravel()
    for ax, col in zip(axes, out.columns):
        ax.plot(out.index, out[col])
        ax.set(title=col, xlabel='Time')
        if ax in [axes[0], axes[4]]:
            ax.set_ylabel('Concentration')
    plt.savefig(os.path.join(here, '../data_out/fig.png'))
    plt.clf()
    timecourse_filename = ('timecourse_'
                           + path_from_here_to_sbml_file
                           .split('/')[-1]
                           .replace('.xml', '.csv'))
    out_target = os.path.join(here, '../data_out/' + timecourse_filename)
    print(f'Writing results to {out_target}...')
    out.to_csv(out_target)
    print('Finished!')
