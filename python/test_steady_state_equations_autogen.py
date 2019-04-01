import arviz
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import pystan
from python_functions import StanModel_cache


TIME_POINTS = np.linspace(0, 100, 100)

ode_metabolites = {
    'DihydroxyacetonePhosphate': 22.2255,
    'PhosphatesInGlycosome': 1.86859,
    'Glyceraldehyde3_phosphate': 0.00976759,
    'Pyruvate': 26.0448,
    'Glucose6_phosphate': 0.506809,
    'PhosphatesCytosol': 35.5604,
    'Glucose': 0.193805,
    'NAD': 0.967427,
    'Fructose6_phosphate': 0.12518,
    'bisphosphoglycerate13': 0.00799218,
    'PGA3PGA2PEP': 9.09737,
    'FructoseBisphosphate16': 4.04497
}
known_reals = {
    'ct_0': 0.97840095394,
    'totalCell': 5.7,
    'glycosome': 0.2446,
    'cytosol': 5.4554,
    'extracellular': 1.0,
    'totVolumePerMgProtein': 5.7,
    'TPIact': 1.0,
    'sumAc': 3.9,
    'sumAg': 6.0,
    'KeqAK': 0.442,
    'Keq_anti': 1.0,
    'sumc4': 45.0,
    'sumc5': 5.0,
    'Keq_PGM': 0.187,
    'Keq_ENO': 6.7,
    'Glycerol': 0.0,
    'PyruvateExternal': 0.0,
    'GlucoseExternal': 5.0
}
kinetic_parameters = {
    'Vm1': 106.2,
    'K1Glc': 2.0,
    'afac': 0.75,
    'Vm2': 625.0,
    'K2ATPg': 0.116,
    'K2GlcI': 0.1,
    'K2Glc6P': 12.0,
    'K2ADPg': 0.126,
    'Vm3': 848.0,
    'K3Glc6P': 0.4,
    'K3Fru6P': 0.12,
    'K4i1Fru16BP': 15.8,
    'Vm4': 780.0,
    'K4ATPg': 0.026,
    'K4Fru6P': 0.82,
    'K4i2Fru16BP': 10.7,
    'Vm5r': 219.555,
    'K5DHAP': 0.015,
    'K5GAP': 0.067,
    'Vm5f': 184.5,
    'K5GAPi': 0.098,
    'Vm6': 842.0,
    'K6GAP': 0.25,
    'K6DHAPg': 1.2,
    'Vm7': 1.0,
    'Vm7f': 1470.0,
    'K7GAP': 0.15,
    'K7NAD': 0.45,
    'Vm7r': 984.9,
    'K7BPGA13': 0.1,
    'K7NADH': 0.02,
    'Vm8': 1.0,
    'Vm8f': 533.0,
    'K8DHAPg': 0.1,
    'K8NADH': 0.01,
    'Vm8r': 149.24,
    'K8Gly3Pg': 2.0,
    'K8NAD': 0.4,
    'Vm9': 368.0,
    'K9Gly3Pc': 1.7,
    'Vm10': 200.0,
    'K10Pyr': 1.96,
    'Vm11': 1.0,
    'Vm11f': 640.0,
    'Vm11r': 18.56,
    'K11ATPg': 0.29,
    'K11PGA3': 1.62,
    'K11ADPg': 0.1,
    'K11BPGA13': 0.05,
    'n12': 2.5,
    'Vm12': 2600.0,
    'K12ADP': 0.114,
    'K13': 50.0,
    'Vm14': 1.0,
    'Vm14r': 33400.0,
    'K14ATPg': 0.19,
    'K14Gly': 0.12,
    'Vm14f': 200.0,
    'K14ADPg': 0.12,
    'K14Gly3Pg': 5.1
}
derived_quantities = [
    'NADH',
    'DihydroxyacetonePhosphate_c',
    'PhosphatesInGlycosome_c',
    'Glyceraldehyde3_phosphate_c',
    'Pyruvate_c',
    'Glucose6_phosphate_c',
    'PhosphatesCytosol_c',
    'Glucose_c',
    'NAD_c',
    'Fructose6_phosphate_c',
    'bisphosphoglycerate13_c',
    'PGA3PGA2PEP_c',
    'FructoseBisphosphate16_c',
    'NADH_c',
    'Glycerol_c',
    'PyruvateExternal_c',
    'GlucoseExternal_c',
    'Vc',
    'Vg',
    'ATPCyt',
    'ATPGly',
    'DHAPCyt',
    'PGA3G',
    'ATPGly_c',
    'PGA3G_c',
    'ATPCyt_c',
    'DHAPCyt_c',
    'DHAPGly',
    'ADPGly',
    'PEPC',
    'ADPCyt',
    'Gy3PC',
    'DHAPGly_c',
    'ADPGly_c',
    'PEPC_c',
    'ADPCyt_c',
    'Gy3PC_c',
    'Gy3PG',
    'Gy3PG_c',
    'Glycerol3_phosphate',
    'Glycerol3_phosphate_c'
]




if __name__ == '__main__':
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
                                     'derived_quantity': derived_quantities},
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
