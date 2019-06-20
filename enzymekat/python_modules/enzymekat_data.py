import pandas as pd
import toml
from typing import Dict, List

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
                 reactions: Dict,
        self.constants = constants
        self.experiments = experiments
        self.reactions = reactions


def from_toml(path):
    t = toml.load(path)
    return EnzymeKatData(
        constants=t['constants'],
        experiments=t['experiments'],
        reactions=t['reactions']
    )
