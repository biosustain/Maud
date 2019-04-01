import arviz
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import pystan
from python_functions import StanModel_cache
from hardcoded_numbers import (EXAMPLE_SPECIES_INPUT,
                               EXAMPLE_KINETIC_PARAMETER_INPUT,
                               EXAMPLE_KNOWN_REALS)
DERIVED_METABOLITES = [
    'PYR', 'eiP', 'hprP', 'NAD', 'AMP', 'BPG', 'eiiaP', 'GLCx', 'eiicbP',
    'MgADP', 'MgATP', 'MgFDP'
]
TIME_POINTS = np.linspace(0, 100, 100)

if __name__ == '__main__':
    model = StanModel_cache(file='test_steady_state_equations.stan')
    data = {
        'N_ode': len(EXAMPLE_SPECIES_INPUT.values()),
        'N_derived': len(DERIVED_METABOLITES),
        'N_known_real': len(EXAMPLE_KNOWN_REALS.values()),
        'P': len(EXAMPLE_KINETIC_PARAMETER_INPUT.values()),
        'T': len(TIME_POINTS) - 1,
        'initial_metabolite_ode': list(EXAMPLE_SPECIES_INPUT.values()),
        'kinetic_parameters': list(EXAMPLE_KINETIC_PARAMETER_INPUT.values()),
        'known_reals': list(EXAMPLE_KNOWN_REALS.values()),
        'ts': TIME_POINTS[1:],
        't0': TIME_POINTS[0]

    }
    fit = model.sampling(data=data, algorithm='Fixed_param', iter=1, chains=1)
    infd = arviz.from_pystan(posterior=fit,
                             coords={'sim_time': TIME_POINTS,
                                     'metabolite_ode': list(EXAMPLE_SPECIES_INPUT.keys()),
                                     'metabolite_derived': DERIVED_METABOLITES},
                             dims={'metabolite_sim_ode': ['sim_time', 'metabolite_ode'],
                                   'metabolite_sim_derived': ['sim_time', 'metabolite_derived']})
    ode_out = infd.posterior['metabolite_sim_ode'].mean(dim=['chain', 'draw']).to_series().unstack()
    derived_out = infd.posterior['metabolite_sim_derived'].mean(dim=['chain', 'draw']).to_series().unstack()
    out = ode_out.join(derived_out)
    out = out.sort_index()
    f, axes = plt.subplots(5, 6, sharex=True, figsize=[15, 15])
    axes = axes.ravel()
    for ax, col in zip(axes, out.columns):
        ax.plot(out.index, out[col])
        ax.set(title=col, xlabel='Time')
        if ax in [axes[0], axes[4]]:
            ax.set_ylabel('Concentration')
    plt.savefig('fig.png')
    plt.clf()
    out.to_csv('ode_species.csv')
    print(out)
