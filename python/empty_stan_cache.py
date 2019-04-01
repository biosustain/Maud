import os

if __name__ == '__main__':
    here = os.path.dirname(os.path.realpath(__file__))
    here_to_there = '../stan/cached_stan_models/'
    there = os.path.join(here, here_to_there)
    targets = os.listdir(there)

    for target in targets:
        if target.endswith(".pkl"):
            os.remove(os.path.join(there, target))
