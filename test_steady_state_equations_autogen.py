import arviz
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import pystan
from convert_mmd_to_stan import get_ode_metabolites, get_known_reals, get_named_quantity, get_derived_quantities
from python_functions import StanModel_cache


MMD_FILE = 'Ecoli glycolisis.mmd'
TIME_POINTS = np.linspace(0, 100, 100)


if __name__ == '__main__':
    ode_metabolites = get_ode_metabolites(MMD_FILE)
    known_reals = get_known_reals(MMD_FILE)
    derived_quantities = get_derived_quantities(MMD_FILE)
    kinetic_parameters = get_named_quantity(MMD_FILE, 'kinetic parameter')
    model = StanModel_cache(file='test_steady_state_autogen.stan')
    data = {
        'N_ode': len(ode_metabolites),
        'N_derived': len(derived_quantities),
        'N_known_real': len(known_reals),
        'P': len(kinetic_parameters.values()),
        'T': len(TIME_POINTS) - 1,
        'initial_metabolite_ode': list(ode_metabolites.values()),
        'kinetic_parameters': list(kinetic_parameters.values()),
        'known_reals': list(known_reals.values()),
        'ts': TIME_POINTS[1:],
        't0': TIME_POINTS[0]

    }
    fit = model.sampling(data=data, algorithm='Fixed_param', iter=1, chains=1)
    infd = arviz.from_pystan(posterior=fit,
                             coords={'sim_time': TIME_POINTS,
                                     'metabolite_ode': list(ode_metabolites.keys()),
                                     'derived_quantity': list(derived_quantities.keys())},
                             dims={'ode_metabolite_sim': ['sim_time', 'metabolite_ode'],
                                   'derived_quantity_sim': ['sim_time', 'derived_quantity']})
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
    plt.savefig('fig.png')
    plt.clf()
    out.to_csv('ode_species.csv')
    print(out)
