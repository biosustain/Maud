import arviz
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import pystan
from python_functions import StanModel_cache
from hardcoded_numbers import (EXAMPLE_SPECIES_INPUT,
                               EXAMPLE_KINETIC_PARAMETER_INPUT,
                               EXAMPLE_KNOWN_REALS)
DERIVED_QUANTITIES = [
    'PYR', 'eiP', 'hprP', 'NAD', 'AMP', 'BPG', 'eiiaP', 'GLCx', 'eiicbP',
    'MgADP', 'MgATP', 'MgFDP'
]
TIME_POINTS = np.linspace(0, 20, 100)

if __name__ == '__main__':
    model = StanModel_cache(file='test_steady_state_equations.stan')
    data = {
        'S': len(EXAMPLE_SPECIES_INPUT.values()),
        'P': len(EXAMPLE_KINETIC_PARAMETER_INPUT.values()),
        'R': len(EXAMPLE_KNOWN_REALS.values()),
        'T': len(TIME_POINTS) - 1,
        'ts': TIME_POINTS[1:],
        't0': TIME_POINTS[0],
        'species': list(EXAMPLE_SPECIES_INPUT.values()),
        'kinetic_parameters': list(EXAMPLE_KINETIC_PARAMETER_INPUT.values()),
        'known_reals': list(EXAMPLE_KNOWN_REALS.values())

    }
    fit = model.sampling(data=data, algorithm='Fixed_param', iter=1, chains=1)
    infd = arviz.from_pystan(posterior=fit,
                             coords={'sim_time': TIME_POINTS,
                                     'species': list(EXAMPLE_SPECIES_INPUT.keys()),
                                     'derived_quantity': DERIVED_QUANTITIES},
                             dims={'species_sim': ['sim_time', 'species'],
                                   'derived_quantities_sim': ['sim_time', 'derived_quantity']})
    species_out = infd.posterior['species_sim'].mean(dim=['chain', 'draw']).to_series().unstack()
    derived_out = infd.posterior['derived_quantities_sim'].mean(dim=['chain', 'draw']).to_series().unstack()
    out = species_out.join(derived_out)
    out = out.sort_index()
    f, axes = plt.subplots(2, 4, sharex=True, figsize=[15, 10])
    axes = axes.ravel()
    for ax, col in zip(axes,
                       ['GLCp', 'ADP', 'ATP', 'P', 'G6P', 'FDP', 'PYR', 'GLCx']):
        ax.plot(out.index, out[col])
        ax.set(title=col, xlabel='Time')
        if ax in [axes[0], axes[4]]:
            ax.set_ylabel('Concentration')
    plt.savefig('fig.png')
    plt.clf()
    out.to_csv('ode_species.csv')
    print(out)
