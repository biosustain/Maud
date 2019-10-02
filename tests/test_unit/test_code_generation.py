import maud.code_generation as code_generation
import maud.io as io
import os

LINEAR_TOML_PATH = os.path.join(os.path.dirname(__file__), '../data/linear.toml')
LINEAR_STAN_CODE_PATH = os.path.join(os.path.dirname(__file__), '../data/linear.stan')


def test_create_inference_program(request):
    mi = io.load_maud_input_from_toml(LINEAR_TOML_PATH)
    correct_stan_code = open(LINEAR_STAN_CODE_PATH, 'r').read()
    generated_stan_code = code_generation.create_stan_program(mi, 'inference')
    assert generated_stan_code == correct_stan_code
    
