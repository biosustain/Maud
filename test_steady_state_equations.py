import pystan
from python_functions import StanModel_cache

SPECIES_INPUT = {}
KINETIC_PARAMETER_INPUT = {}
KNOWN_REALS_INPUT = {}
KNOWN_INTS_INPUT = {}


if __name__ == '__main__':
    model = StanModel_cache(file='test_steady_state_equations.stan')
    data = {
        'S': len(SPECIES_INPUT.values()),
        'P': len(KINETIC_PARAMETER_INPUT.values()),
        'R': len(KNOWN_REALS_INPUT.values()),
        'I': len(KNOWN_INTS_INPUT.values()),
        'species': SPECIES_INPUT.values(),
        'kinetic_parameters': KINETIC_PARAMETER_INPUT.values(),
        'known_reals': KNOWN_REALS_INPUT.values(),
        'known_ints': KNOWN_INTS_INPUT.values()
    }
    fit = model.sampling(data=data, algorithm='Fixed_param')
    for species, dsdt in zip(SPECIES_INPUT.keys(), fit['dsdt'].reshape(-1)):
        print(f'{species}: {dsdt.round(3)}')
