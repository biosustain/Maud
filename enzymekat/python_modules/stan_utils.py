import numpy as np
import os
from multiprocessing import Pool


RELATIVE_PATH_CMDSTAN = '../../cmdstan'


def compile_stan_model_with_cmdstan(path_from_cmdstan_home_to_program):
    here = os.path.dirname(os.path.abspath(__file__))
    cmdstan_home = os.path.join(here, RELATIVE_PATH_CMDSTAN)
    command = f"cd {cmdstan_home} && make {path_from_cmdstan_home_to_program}"
    os.system(command)


def run_compiled_cmdstan_model(program_path,
                               input_path,
                               output_path,
                               method_config='sample',
                               init_config = '',
                               random_config=None,
                               refresh_config='',
                               chains=1):
    """Get a compiled cmdstan program to produce a csv of samples
    
    See the cmdstan manual for more about the commands and how to customize
    them.
    
    The program will write a csv to `output_path` that can be read using
    arviz's `from_pystan` function.

    :param program_path: path to the compiled program :param input_path: path
    to a stan Rdump file of input data 

    :param output_path: path to the output csv 

    :param method_config: customize the command at the method level 

    :param init_config: optionally give the model a path to an Rdump of initial
    parameter values 

    :param random_config: optionally fix the random seed

    :param refresh_config: set the refresh interval for progress reporting

    :param chains: how many parallel chains do you want to run (probably no
    more than there are cores on your computer)
    
    :return: None

    """
    here = os.path.dirname(os.path.abspath(__file__))
    rng = np.random.RandomState()
    program_directory, program = os.path.split(program_path)
    path_template = output_path.replace('.csv', '{prefix}.csv')
    if random_config is None:
        random_config = 'random seed={random}'
    command_template = f"""
        cd {program_directory}

        ./{program} \
        {method_config} \
        data file={input_path} \
        {init_config} \
        {random_config} \
        output file={path_template} \
        {refresh_config}
        """
    commands = []
    for c in range(chains):
        pref = str(c)
        rand = rng.randint(1, 99999 + 1)
        command = command_template.format(prefix=pref, random=rand)
        commands.append(command)

    pool = Pool(processes=chains)
    pool.map(os.system, commands)
    replace = path_template.format(prefix='*')
    cat_command = f"cat {replace} > {output_path}"
    os.system(cat_command)
    
