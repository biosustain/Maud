"""
    Script for running a timecourse based on an sbml file and a stan model
    generated from it.

    Running the timecourse model requires integrating a stiff ODE system, which
    isn't possible with pystan at the moment due to a bug (see
    https://github.com/stan-dev/pystan/pull/518). As a workaround, this script
    uses cmdstan instead, which works but requires some configuration.

    To run this script you will need to follow these steps:

    1. Download (or `git clone`) a recent cmdstan 
       (https://github.com/stan-dev/cmdstan/releases)
    2. Run the following commands from the root of the cmdstan directory
       - `make stan-update`
       - `make clean-all`
       - `make build`
    3. Create a directory for models e.g. `mkdir models` 
    4. Create a directory for the specific timecourse model you want to 
       run, e.g. `mkdir models/t_brucei`
    5. Configure the upper case variables at the top of this script so that they 
       point at the things they should do.

    To build a pdf cmdstan manual, run `make manual`. You can also just run `make`
    to get a broad overview of how cmdstan works.

"""
import arviz
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import pystan
import os
import stan_utils
from sbml_functions import read_sbml_file

# These need to be configured!
PATH_FROM_HERE_TO_CMDSTAN_HOME = '../cmdstan'
PATH_FROM_HERE_TO_STAN_MODEL = '../stan/autogen/t_brucei_timecourse.stan'
PATH_FROM_HERE_TO_SBML_FILE = '../data_in/t_brucei.xml'
PATH_FROM_CMDSTAN_HOME_TO_KINETICS_FILE = '../stan/autogen/t_brucei.stan'
PATH_FROM_HERE_TO_TIMECOURSE_TEMPLATE = '../stan/timecourse_model_template.stan'
PATH_FROM_HERE_TO_CMDSTAN_INPUT_DATA_FILE = '../data_out/timecourse_input.Rdump'
PATH_FROM_HERE_TO_CMDSTAN_OUTPUT_DATA_FILE = '../data_out/t_brucei_timecourse_cmdstan.csv'
PATH_FROM_HERE_TO_TIMECOURSE_CSV = '../data_out/t_brucei_timecourse.csv'
TIME_POINTS = np.linspace(0, 0.5, 101)

if __name__ == '__main__':
    # get some paths
    here = os.path.dirname(os.path.abspath(__file__))
    cmdstan_directory = os.path.join(here, PATH_FROM_HERE_TO_CMDSTAN_HOME)
    template_path = os.path.join(here, PATH_FROM_HERE_TO_TIMECOURSE_TEMPLATE)
    stan_model_path = os.path.join(here, PATH_FROM_HERE_TO_STAN_MODEL)
    cmdstan_output_data_file = os.path.join(here, PATH_FROM_HERE_TO_CMDSTAN_OUTPUT_DATA_FILE)
    cmdstan_input_data_file = os.path.join(here, PATH_FROM_HERE_TO_CMDSTAN_INPUT_DATA_FILE)
    timecourse_csv_path = os.path.join(here, PATH_FROM_HERE_TO_TIMECOURSE_CSV)

    # parse sbml file
    m = read_sbml_file(os.path.join(here, PATH_FROM_HERE_TO_SBML_FILE))

    # save input data for model
    data = {
        'N_ode': len(m.ode_metabolites.values()),
        'N_derived': len(m.assignment_expressions.values()),
        'N_known_real': len(m.known_reals.values()),
        'P': len(m.kinetic_parameters.values()),
        'T': len(TIME_POINTS) - 1,
        'initial_metabolite_ode': list(m.ode_metabolites.values()),
        'kinetic_parameters': list(m.kinetic_parameters.values()),
        'known_reals': list(m.known_reals.values()),
        'ts': TIME_POINTS[1:],
        't0': TIME_POINTS[0],
        'rel_tol': 1e-8,
        'abs_tol': 1e-8,
        'max_num_steps': int(1e8)
    }
    pystan.misc.stan_rdump(data, cmdstan_input_data_file)

    # write timecourse model based on template
    with open(template_path, 'r') as f:
        model_code = f.read().replace(
            'REPLACE_THIS_WORD', PATH_FROM_CMDSTAN_HOME_TO_KINETICS_FILE
        )
        f.close()
        with open(stan_model_path, 'w') as f:
            f.write(model_code)
            f.close()

    # compile model with cmdstan
    if not os.path.isfile(stan_model_path.replace('.stan', '.hpp')):
        compile_path =  os.path.relpath(stan_model_path.replace('.stan', ''),
                                        start=cmdstan_directory)
        stan_utils.compile_stan_model_with_cmdstan(compile_path)

    method_config = 'sample algorithm=fixed_param num_warmup=0 num_samples=1'
    stan_utils.run_compiled_cmdstan_model(
        stan_model_path.replace('.stan', ''),
        cmdstan_input_data_file,
        cmdstan_output_data_file,
        method_config=method_config
    )
    infd = arviz.from_cmdstan([cmdstan_output_data_file],
                              coords={'sim_time': TIME_POINTS,
                                      'metabolite_ode': list(m.ode_metabolites.keys()),
                                      'metabolite_derived': list(m.assignment_expressions.keys())},
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
    print(f'Writing results to {timecourse_csv_path}...')
    out.to_csv(timecourse_csv_path)
    print('Finished!')
