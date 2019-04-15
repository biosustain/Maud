import os

PATH_FROM_HERE_TO_CMDSTAN_HOME = '../cmdstan'


def compile_stan_model_with_cmdstan(path_from_cmdstan_home_to_program):
    here = os.path.dirname(os.path.abspath(__file__))
    cmdstan_home = os.path.join(here, PATH_FROM_HERE_TO_CMDSTAN_HOME)
    command = f"cd {cmdstan_home} && make {path_from_cmdstan_home_to_program}"
    os.system(command)


def run_compiled_cmdstan_model(program_path,
                               input_path,
                               output_path,
                               method_config='sample',
                               init_config = '',
                               random_config='',
                               refresh_config=''):
    here = os.path.dirname(os.path.abspath(__file__))
    
    program_directory, program = os.path.split(program_path)
    print(program)

    command = f"""
    cd {program_directory}

    ./{program} \
    {method_config} \
    data file={input_path} \
    {init_config} \
    {random_config} \
    output file={output_path} \
    {refresh_config}
    """
    print(command)
    os.system(command)
