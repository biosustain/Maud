import pandas as pd
import toml
from typing import Dict, List
from python_modules.stan_utils import one_index_factorize


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
        self.stoichiometry = self.get_stoichoimetry()
        self.measurements = self.get_measurements()
        self.known_reals = self.get_known_reals()
        self.parameters = self.get_parameters()
        self.metabolites = self.get_metabolites()

    def get_stoichoimetry(self):
        out = (
            pd.DataFrame.from_records(self.reactions)
            .set_index('name')
            ['stoichiometry']
            .pipe(expand_series_of_dicts)
            .fillna(0)
            .astype(int)
        )
        out.index.name = 'reaction'
        out.columns.name = 'metabolite'
        return out

    def get_measurements(self):
        out = pd.io.json.json_normalize(
            self.experiments,
            'measurements',
            meta=['label'],
            meta_prefix='experiment_'
        ).rename(columns={'experiment_label': 'experiment'})
        out['metabolite'] = out['label'].where(out['type'] == 'metabolite')
        out['reaction'] = out['label'].where(out['type'] == 'flux')
        out['metabolite_code'], _ = one_index_factorize(out['metabolite'])
        out['reaction_code'], _ = one_index_factorize(out['reaction'])
        out['experiment_code'], _ = one_index_factorize(out['experiment'])
        return out

    def get_known_reals(self):
        out = {}
        for e in self.experiments:
            out[e['label']] = {**e['conditions'], **self.constants}
        return pd.DataFrame.from_dict(out).assign(stan_code=lambda df: range(1, len(df) + 1))

    def get_parameters(self):
        return (
            pd.io.json.json_normalize(self.reactions,
                                      'parameters',
                                      meta='name')
            .assign(stan_code=lambda df: range(1, len(df) + 1))
            .rename(columns={'name': 'reaction',
                             'label': 'parameter'})
        )

    def get_metabolites(self):
        return (
            pd.DataFrame({'name': self.stoichiometry.columns})
            .assign(stan_code=lambda df: range(1, len(df) + 1))
        )
            

def from_toml(path):
    t = toml.load(path)
    return EnzymeKatData(
        constants=t['constants'],
        experiments=t['experiments'],
        reactions=t['reactions']
    )
