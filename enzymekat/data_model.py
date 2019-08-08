import pandas as pd
from pandas.io.json import json_normalize
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
                 experiments: List[Dict],
                 reactions: Dict):
        self.constants = constants
        self.experiments = experiments
        self.reactions = reactions
        self.stoichiometry = self.get_stoichiometry()
        self.stan_codes = self.get_stan_codes()
        self.metabolites = self.get_metabolites()
        self.flux_measurements, self.concentration_measurements = self.get_measurements()
        self.known_reals = self.get_known_reals()
        self.kinetic_parameters = self.get_kinetic_parameters()
        self.unbalanced_metabolite_priors = self.get_unbalanced_metabolite_priors()

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

    def get_stan_codes(self):

        def map_to_codes(l: List[str]):
            return dict(zip(l, range(1, len(l) + 1)))

        S = self.get_stoichiometry()
        experiments = self.experiments
        reactions = self.reactions
        reaction_names = list(S.index.unique())
        metabolite_names = list(S.columns.unique())
        experiment_names = [e['label'] for e in experiments]
        kinetic_parameter_names = []
        for r in reactions:
            for p in r['parameters']:
                name = r['name'] + '_' + p['label']
                if 'metabolite' in p.keys():
                    name = name + '_' + p['metabolite']
                kinetic_parameter_names.append(name)
        return {
            'reaction': map_to_codes(reaction_names),
            'metabolite': map_to_codes(metabolite_names),
            'experiment': map_to_codes(experiment_names),
            'kinetic_parameter': map_to_codes(kinetic_parameter_names)
        }

    def get_unbalanced_metabolite_names(self):
        S = self.get_stoichiometry()
        produced_but_not_consumed = S.gt(0).any(axis=0) & S.ge(0).all(axis=0)
        consumed_but_not_produced = S.lt(0).any(axis=0) & S.le(0).all(axis=0)
        is_unbalanced = produced_but_not_consumed | consumed_but_not_produced
        return S.columns[is_unbalanced]

    def get_unbalanced_metabolite_priors(self):
        stan_codes = self.get_stan_codes()
        experiments = self.experiments
        metabolite_codes = pd.Series(stan_codes['metabolite'], name='metabolite_code')
        experiment_codes = pd.Series(stan_codes['experiment'], name='experiment_code')
        priors = (
            pd.io.json.json_normalize(experiments,
                                      record_path='unbalanced_metabolite_priors',
                                      meta='label',
                                      meta_prefix='experiment_')
            .rename(columns={'label': 'metabolite'})
        )
        # correct dtypes - they are object but should be str
        object_columns = ['metabolite', 'experiment_label']
        priors[object_columns] = priors[object_columns].astype(str)
        return (
            priors
            .join(metabolite_codes, on='metabolite')
            .join(experiment_codes, on='experiment_label')
        )

    def get_metabolites(self):
        S = self.get_stoichiometry()
        metabolite_names = S.columns
        stan_codes = pd.Series(self.get_stan_codes()['metabolite'], name='stan_code')
        unbalanced_metabolite_names = self.get_unbalanced_metabolite_names()
        is_unbalanced = [m in unbalanced_metabolite_names for m in metabolite_names]
        return pd.DataFrame({
            'name': metabolite_names,
            'is_unbalanced': is_unbalanced
        }).join(stan_codes, on='name')

    def get_measurements(self):
        stan_codes = self.get_stan_codes()
        experiments = self.experiments
        metabolite_codes = stan_codes['metabolite']
        reaction_codes = stan_codes['reaction']
        experiment_codes = stan_codes['experiment']
        measurements = pd.io.json.json_normalize(
            experiments,
            'measurements',
            meta=['label'],
            meta_prefix='experiment_'
        )
        flux_measurements, conc_measurements = (
            measurements
            .query(f"type == '{t}'")
            .dropna(how='all', axis=1)
            .copy()
            .rename(columns={'label': t + '_label'})
            for t in ['flux', 'metabolite']
        )
        flux_measurements['reaction_code'] = flux_measurements['flux_label'].map(reaction_codes.get)
        flux_measurements['experiment_code'] = flux_measurements['experiment_label'].map(experiment_codes.get)
        conc_measurements['metabolite_code'] = conc_measurements['metabolite_label'].map(metabolite_codes.get)
        conc_measurements['experiment_code'] = conc_measurements['experiment_label'].map(experiment_codes.get)
        return flux_measurements, conc_measurements

    def get_known_reals(self):
        out = {}
        for e in self.experiments:
            conditions = {k: v for k, v in e.items() if type(v) in [float, int]}
            out[e['label']] = {**conditions, **self.constants}
        return pd.DataFrame.from_dict(out).assign(stan_code=lambda df: range(1, len(df) + 1))

    def get_kinetic_parameters(self):
        reactions = self.reactions
        stan_codes = pd.Series(self.get_stan_codes()['kinetic_parameter'], name='stan_code')
        for r in reactions:  # json_normalize doesn't work with missing record path
            if 'parameters' not in r.keys():
                r['parameters'] = []
        pars = (
            json_normalize(reactions, record_path='parameters', meta='name')
            .rename(columns={'name': 'reaction', 'label': 'parameter'})
            .assign(is_allosteric=lambda df: df['type'].eq('allosteric').astype(int))
        )
        label_cols = (
            ['reaction', 'parameter', 'metabolite']
            if 'metabolite' in pars.columns
            else ['reaction', 'parameter']
        )
        pars['label'] = pars[label_cols].apply(lambda row: '_'.join(row.dropna().astype(str)), axis=1)
        pars = pars.join(stan_codes, on='label')
        return pars


def from_toml(path):
    t = toml.load(path)
    return EnzymeKatData(
        constants=t['constants'],
        experiments=t['experiment'],
        reactions=t['reaction']
    )


def join_string_values(df):
    return df.apply(lambda x: '_'.join(x.dropna().astype(str)), axis=1)
