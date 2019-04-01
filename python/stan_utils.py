import pystan
import pickle
from hashlib import md5


def StanModel_cache(file, model_name=None, **kwargs):
    with open(file, 'r') as f:
        model_code = f.read()
        f.close()
    code_hash = md5(model_code.encode('ascii')).hexdigest()
    if model_name is None:
        cache_fn = 'data/cached_stan_models/{}.pkl'.format(code_hash)
    else:
        cache_fn = 'data/cached_stan_models/{}-{}.pkl'.format(model_name, code_hash)
    try:
        sm = pickle.load(open(cache_fn, 'rb'))
    except:
        sm = pystan.StanModel(model_code=model_code)
        with open(cache_fn, 'wb') as f:
            pickle.dump(sm, f)
    else:
        print("Using cached StanModel")
    return sm
