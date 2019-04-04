import os
import pystan
import pickle
from hashlib import md5


def StanModel_cache(model_input, model_name=None, **kwargs):
    here = os.path.dirname(__file__)
    with open(model_input, 'r') as f:
        model_code = f.read()
        f.close()
    code_hash = md5(model_code.encode('ascii')).hexdigest()
    if model_name is None:
        cache_fn = os.path.join(here, '../stan/cached_stan_models/{}.pkl'.format(code_hash))
    else:
        cache_fn = os.path.join(here, '../stan/cached_stan_models/{}-{}.pkl'.format(model_name, code_hash))
    try:
        sm = pickle.load(open(cache_fn, 'rb'))
    except:
        sm = pystan.StanModel(model_code=model_code)
        with open(cache_fn, 'wb') as f:
            pickle.dump(sm, f)
    else:
        print("Using cached StanModel")
    return sm


def run_cmdstan_model(cmdstan_directory,
                      path_from_cmdstan_directory_to_program,
                      cmdstan_input_data_file,
                      cmdstan_output_data_file):

    here = os.path.dirname(__file__)
    os.environ["CMDSTAN_DIRECTORY"] = cmdstan_directory
    os.environ["PATH_FROM_CMDSTAN_DIRECTORY_TO_PROGRAM"] = path_from_cmdstan_directory_to_program
    os.environ["CMDSTAN_INPUT_DATA_FILE"] = cmdstan_input_data_file
    os.environ["CMDSTAN_OUTPUT_DATA_FILE"] = cmdstan_output_data_file

    if not os.path.isfile(os.path.join(cmdstan_directory, path_from_cmdstan_directory_to_program)):
        os.system("bash " + os.path.join(here, "../bash/compile_timecourse_model.sh"))

    os.system("bash " + os.path.join(here, "../bash/run_timecourse_model.sh"))
