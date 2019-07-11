import pandas as pd
import toml
from typing import Dict, List


def expand_series_of_dicts(s):
    """get a dataframe from a series whose values are dictionaries"""
    return pd.DataFrame.from_records(s.values, index=s.index)



class EnzymeKatData():
    """Object with all the data that is needed to fit an enzymeKat model, namely:

      - Information about the experimental setup (measurement accuracy, time to
        steady state)
      - Information about kinetic parameters
      - Information about ode metabolites

    """
    def __init__(self,
                 constants: Dict,
                 experiments: Dict,
                 reactions: Dict):
        self.constants = constants
        self.experiments = experiments
        self.reactions = reactions
        self.stoichiometry = self.get_stoichiometry()
        self.flux_measurements, self.concentration_measurements = self.get_measurements()
        self.known_reals = self.get_known_reals()
        self.kinetic_parameters, self.thermodynamic_parameters = self.get_parameters()
        self.stan_codes = self.get_stan_codes()

    def get_stoichiometry(self):
        S = (
            pd.DataFrame.from_records(self.reactions)
            .set_index('name')
            ['stoichiometry']
            .pipe(expand_series_of_dicts)
            .fillna(0)
            .astype(int)
        )
        S.index.name = 'reaction'
        S.columns.name = 'metabolite'
        return S

    def get_measurements(self):
        stan_codes = self.get_stan_codes()
        measurements = pd.io.json.json_normalize(
            self.experiments,
            'measurements',
            meta=['label'],
            meta_prefix='experiment_'
        )
        measurements['reaction_code'] = measurements['label'].map(stan_codes['reaction'].get)
        measurements['metabolite_code'] = measurements['label'].map(stan_codes['metabolite'].get)
        measurements['experiment_code'] = measurements['experiment_label'].map(stan_codes['experiment'].get)
        flux_measurements, conc_measurements = (
            measurements
            .query(f"type == '{t}'")
            .dropna(how='all', axis=1)
            .copy()
            .rename(columns={'label': t + '_label'})
            for t in ['flux', 'metabolite']
        )
        flux_measurements['reaction_code'] = flux_measurements['reaction_code'].astype(int)
        conc_measurements['metabolite_code'] = conc_measurements['metabolite_code'].astype(int)
        return flux_measurements, conc_measurements

    def get_stan_codes(self):
        S = self.get_stoichiometry()
        known_reals = self.get_known_reals()
        kinetic_parameters, thermodynamic_parameters = self.get_parameters()
        experiment_labels = [e['label'] for e in self.experiments]
        reaction_codes = dict(zip(S.index, range(1, len(S.index) + 1)))
        metabolite_codes = dict(zip(S.columns, range(1, len(S.columns) + 1)))
        experiment_codes = dict(zip(experiment_labels, range(1, len(experiment_labels) + 1)))
        known_real_codes = dict(zip(known_reals.index, range(1, len(known_reals.index) + 1)))
        kinetic_parameter_codes = dict(zip(kinetic_parameters['label'],
                                           range(1, len(kinetic_parameters) + 1)))
        thermodynamic_parameter_codes = dict(zip(thermodynamic_parameters['label'],
                                                 range(1, len(kinetic_parameters) + 1)))
        return {'reaction': reaction_codes,
                'metabolite': metabolite_codes,
                'experiment': experiment_codes,
                'known_real': known_real_codes,
                'kinetic_parameter': kinetic_parameter_codes,
                'thermodynamic_parameter': thermodynamic_parameter_codes}

    def get_known_reals(self):
        out = {}
        for e in self.experiments:
            out[e['label']] = {**e['conditions'], **self.constants}
        return pd.DataFrame.from_dict(out).assign(stan_code=lambda df: range(1, len(df) + 1))

    def get_parameters(self):
        parameters = (
            pd.io.json.json_normalize(self.reactions,
                                      'parameters',
                                      meta='name')
            .rename(columns={'name': 'reaction', 'label': 'parameter'})
            .assign(label=lambda df: df['reaction'].str.cat(df['parameter'], sep='_'))
        )
        kinetic_parameters = parameters.query("type == 'kinetic'").copy()
        thermodynamic_parameters = parameters.query("type == 'thermodynamic'").copy()
        return kinetic_parameters, thermodynamic_parameters


def from_toml(path):
    t = toml.load(path)
    return EnzymeKatData(
        constants=t['constants'],
        experiments=t['experiments'],
        reactions=t['reactions']
    )
