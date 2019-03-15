import arviz
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import pystan
from python_functions import StanModel_cache

SPECIES_INPUT = {
    'ATP': 2.5722,
    'PEP': 0.997038,
    'P': 9.76395,
    'GAP': 0.117183,
    'F6P': 0.261766,
    'DAP': 0.437094,
    'eiia': 0.0142019,
    'GLCp': 0.00403337,
    'PGA2': 0.378297,
    'ei': 0.000334013,
    'PGA3': 0.696274,
    'eiicb': 4.7191e-05,
    'FDP': 0.281808,
    'hpr': 0.000191212,
    'ADP': 0.598315,
    'G6P': 0.86113,
    'NADH': 0.158456,
}
KINETIC_PARAMETER_INPUT = {
    'Keq': 0.36,
    'KmF6P': 0.147,
    'KmG6P': 0.28,
    'KmPEP': 1.999,
    'Vmax': 2.32456,
    'KmPGN': 0.515958,
    'KefrADP': 0.0735264,
    'KefrPEP': 19.98,
    'KeftADP': 9.009,
    'KeftPEP': 0.26026,
    'Keq_1': 1998,
    'KirADP': 54.945,
    'KirATP': 2.4975e-05,
    'KirF6P': 1.84615,
    'KirFDP': 0.045954,
    'KitADP': 80.08,
    'KitATP': 0.014014,
    'KitF6P': 0.00856856,
    'KitFDP': 50.5505,
    'KmrADP': 0.690009,
    'KmrATPMg': 8.12187e-05,
    'KmrF6P': 2.05205e-05,
    'KmrFDP': 10.01,
    'KmtADP': 2.002,
    'KmtATPMg': 3.34334,
    'KmtF6P': 32.967,
    'KmtFDP': 9.99,
    'L0': 14.0851,
    'Vmax_1': 0.185253,
    'Wr': 0.0237041,
    'Wt': 0.146735,
    'n': 4,
    'Keq_2': 0.18981,
    'KmDAP': 0.13001,
    'KmFDP': 0.12012,
    'KmGAP': 0.13001,
    'KmPEP_1': 0.5,
    'Vmax_2': 21.6978,
    'Keq_3': 0.270203,
    'KmDAP_1': 0.01,
    'KmGAP_1': 1.89301,
    'Vmax_3': 24.1843,
    'Keq_4': 20,
    'KmBPG': 0.2,
    'KmGAP_2': 2.47265,
    'KmNAD': 0.0110454,
    'KmNADH': 3.69797,
    'KmP': 0.017,
    'Vmax_4': 8.66573,
    'Keq_5': 99.9925,
    'KmADPMg': 0.085416,
    'KmATPMg': 3.47737,
    'KmBPG_1': 0.0113296,
    'KmPGA3': 2.45722,
    'Vmax_5': 16.1089,
    'Keq_6': 0.565818,
    'KmPGA2': 1.9153,
    'KmPGA3_1': 0.115,
    'Vmax_6': 10.9934,
    'Keq_7': 3,
    'KmPEP_2': 0.1,
    'KmPGA2_1': 0.1,
    'Vmax_7': 11.7189,
    'KefrFDP': 0.449149,
    'KefrG6P': 0.158746,
    'KefrGL6P': 0.150482,
    'KefrR5P': 9.33254,
    'KefrRU5P': 1.53591,
    'KefrS7P': 0.0785955,
    'KefrX5P': 0.677374,
    'KeftATP': 3.69117,
    'KeftSUCCOA': 8.26406,
    'KirADP_1': 0.517585,
    'KirATP_1': 96.0333,
    'KirPEP': 0.181056,
    'KirPYR': 15.1403,
    'KirPyrATP': 230.781,
    'KitADP_1': 0.224911,
    'KitATP_1': 0.039564,
    'KitPEP': 0.465672,
    'KitPYR': 0.2499,
    'KitPyrATP': 11.3691,
    'KmrADPMg': 0.326144,
    'KmrPEP': 5.56368e-07,
    'KmtADPMg': 0.054678,
    'KmtPEP': 0.11475,
    'L0_1': 50.4818,
    'Vmax_8': 0.74716,
    'n_1': 4,
    'KirAMP': 0.00122122,
    'KirAMPFDP': 0.256256,
    'KirF6P_1': 1.12112,
    'KirF6PMg': 0.384615,
    'KirFDP_1': 1.35327,
    'KirFDPMg': 0.75924,
    'KirFDPMgMg': 0.356356,
    'KirP': 3.16316,
    'KirPF6P': 6.60538,
    'KirPF6PMg': 48.4484,
    'KirPMg': 0.856,
    'KitAMP': 0.000255,
    'KitAMPFDP': 690,
    'KitF6P_1': 0.304,
    'KitF6PMg': 315,
    'KitFDP_1': 0.043101,
    'KitFDPMg': 0.00642,
    'KitFDPMgMg': 100,
    'KitP': 0.642,
    'KitPF6P': 0.00689,
    'KitPF6PMg': 16.5,
    'KitPMg': 539,
    'KmrFDP_1': 0.0636141,
    'KmrMg': 0.039039,
    'KmtFDP_1': 1e-05,
    'KmtMg': 55.055,
    'L0_2': 0.000815,
    'Vmax_9': 0.215583,
    'n_2': 4,
    'KdAMP': 1480,
    'KdATPMgPPS': 0.0549,
    'KdMg': 36.9,
    'KdP': 346,
    'KdPEP': 95.7,
    'KdPYR': 2740,
    'KefADP': 0.0283,
    'KefAKG': 0.274,
    'KefATP': 0.000628,
    'KefOAA': 0.796,
    'Keq_8': 200000,
    'KmAMP': 0.000384,
    'KmATPMg_1': 0.0549,
    'KmP_1': 84.4,
    'KmPEP_3': 20.7,
    'KmPYR': 0.229,
    'Vmax_10': 0.0163772,
    'W': 10,
    'alpha': 38900,
    'KmPEP_4': 0.6,
    'KmPYR_1': 1,
    'kF': 12000,
    'kR': 8000,
    'k1': 200000,
    'k2': 8000,
    'k1_1': 61000,
    'k2_1': 47000,
    'k1_2': 11000,
    'k2_2': 4000,
    'KmG6P_1': 2125.91,
    'KmGLC': 0.02,
    'kF_1': 4000,
    'kR_1': 1e-05,
    'Vmax_11': 1.30166,
    'Keq_9': 3.63369,
    'Vmax_12': 100,
    'Km': 10
}
KNOWN_REALS_INPUT = {
    'PYR': -7.7521005040965605,
    'eiP': 0.0067151757729566096,
    'hprP': 0.0054420236094104221,
    'NAD': 1.5700000006962138,
    'AMP': 3.356768777029624,
    'BPG': 9.9675557170810123,
    'eiiaP': 0.4947703996764079,
    'GLCx': 2.2651042736100142,
    'eiicbP': 0.00038981236144815099,

    'AKG': 0.59787,
    'GL6P': 0.00326165,
    'OAA': 0.12784,
    'PGN': 0.1316,
    'R5P': 0.106842,
    'RU5P': 0.341827,
    'S7P': 0.141985,
    'SUCCOA': 0.0410878,
    'X5P': 0.506018,
 
}
KNOWN_INTS_INPUT = {}
TIME_POINTS = np.linspace(0, 20, 100)

if __name__ == '__main__':
    model = StanModel_cache(file='test_steady_state_equations.stan')
    data = {
        'S': len(SPECIES_INPUT.values()),
        'P': len(KINETIC_PARAMETER_INPUT.values()),
        'R': len(KNOWN_REALS_INPUT.values()),
        'T': len(TIME_POINTS) - 1,
        'ts': TIME_POINTS[1:],
        't0': TIME_POINTS[0],
        'species': list(SPECIES_INPUT.values()),
        'kinetic_parameters': list(KINETIC_PARAMETER_INPUT.values()),
        'known_reals': list(KNOWN_REALS_INPUT.values())

    }
    fit = model.sampling(data=data, algorithm='Fixed_param', iter=1, chains=1)
    infd = arviz.from_pystan(posterior=fit,
                             coords={'sim_time': TIME_POINTS[1:],
                                     'species': list(SPECIES_INPUT.keys())},
                             dims={'species_sim': ['sim_time', 'species']})
    out = infd.posterior['species_sim'].mean(dim=['chain', 'draw']).to_series().unstack()
    out.loc[0] = pd.Series(SPECIES_INPUT)
    out = out.sort_index()
    f, axes = plt.subplots(2, 3, sharex=True)
    axes = axes.ravel()
    for ax, col in zip(axes,
                       ['GLCp', 'ADP', 'ATP', 'P', 'G6P', 'FDP']):
        ax.plot(out.index, out[col])
        ax.set(title=col, xlabel='Time', ylabel='Concentration')
    plt.savefig('fig.png')
    plt.clf()
    out.to_csv('ode_species.csv')
    print(out)