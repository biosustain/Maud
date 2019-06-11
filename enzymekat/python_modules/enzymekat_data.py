import pandas as pd
import toml

class EnzymeKatData():
    """Object with all the data that is needed to fit an enzymeKat model, namely:

      - Information about the experimental setup (measurement accuracy, time to
        steady state)
      - Information about kinetic parameters
      - Information about ode metabolites

    """
    def __init__(self,
                 experiment_info: pd.Series,
                 ode_metabolites: pd.Series,
                 known_reals: pd.Series,
                 reactions: pd.DataFrame,
                 thermodynamic_parameters: pd.DataFrame,
                 kinetic_parameters: pd.DataFrame,
                 enzyme_concentrations: pd.DataFrame):
        self.experiment_info = experiment_info
        self.ode_metabolites = ode_metabolites
        self.known_reals = known_reals
        self.reactions = reactions
        self.thermodynamic_parameters = thermodynamic_parameters
        self.kinetic_parameters = kinetic_parameters
        self.enzyme_concentrations = enzyme_concentrations


def from_toml(path):
    t = toml.load(path)
    return EnzymeKatData(
        experiment_info=pd.Series(t['experiment_info']),
        ode_metabolites=pd.DataFrame.from_records(t['ode_metabolites']),
        known_reals=pd.Series(t['known_reals']),
        reactions=pd.DataFrame.from_records(t['reactions']),
        thermodynamic_parameters=pd.DataFrame.from_records(t['thermodynamic_parameters']),
        kinetic_parameters=pd.DataFrame.from_records(t['kinetic_parameters']),
        enzyme_concentrations=pd.DataFrame.from_records(t['enzyme_concentrations'])
    )
