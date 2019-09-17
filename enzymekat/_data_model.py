"""Definitions of Enzymekat-specific objects"""

from collections import defaultdict
from typing import Dict, List

class Compartment:
    def __init__(self, id : str,
                 name: str = None,
                 volume: float = 1.0):
        """
        Constructor for compartment objects.
​
        :param id: compartment id, use a BiGG id if possible.
        :param name: compartment name.
        :param volume: compartment volume.
        """
        self.id = id
        self.name = name if name is not None else id
        self.volume = volume


class Metabolite:
    def __init__(self, id: str,
                 name: str = None,
                 balanced: bool = None,
                 compartment: Compartment = None):
        """
        Constructor for metabolite objects.
​
        :param id: metabolite id, use a BiGG id if possible.
        :param name: metabolite name.
        :param balanced: Doe this metabolite have an unchanging concentration at steady state?
        :param compartment: compartment for the metabolite.
        """
        self.id = id
        self.name = name if name is not None else id
        self.balanced = balanced
        self.compartment = compartment


class Modifier:
    def __init__(self,
                 metabolite: Metabolite,
                 modifier_type: str = None):
        """
        Constructor for modifier objects.
        
        :param met: the metabolite that is the modifier 
        :param allosteric: whether or not the modifier is allosteric
        :param modifier_type: what is the modifier type: 'activator', 'inhibitor', 'competitive inhibitor',
        'uncompetitive inhibitor', or 'noncompetitive inhibitor' 
        """
        self.metabolite = metabolite
        self.allosteric = modifier_type in ['inhibitor']
        self.modifier_type = modifier_type


class Parameter:
    def __init__(self,
                 id: str,
                 reaction_id: str,
                 metabolite_id: str = None):
        """
        Constructor for parameter object.
        
        :param id: parameter id
        :param reaction_id: id of the reaction associated with the parameter
        :param metabolite_id: id of the metabolite associated with the parameter if any

        """
        self.id = id
        self.reaction_id = reaction_id
        self.metabolite_id = metabolite_id


class Reaction:
    def __init__(self, id: str,
                 name: str = None,
                 reversible: bool = True,
                 is_exchange: bool = None,
                 stoichiometry: Dict = defaultdict(),
                 modifiers: Dict[str, Modifier] = defaultdict(),
                 parameters: Dict[str, Parameter] = defaultdict(),
                 rate_law: str = None,
                 enzymes: set = None):
        """
        Constructor for the reaction object.

        :param id: reaction id, use a BiGG id if possible.
        :param name: reaction name.
        :param reversible: whether or not reaction is reversible.
        :param is_exchange: whether or not reaction is an exchange reaction.
        :param stoichiometry: reaction stoichiometry, e.g. for the reaction: 1.5 f6p <-> fdp we have {'f6p'; -1.5, 'fdp': 1}
        :param modifiers: reaction modifiers, given as {'modifier_id': modifier_object}
        :param parameters: reaction parameters, give as {'parameter_id', parameter_object}
        :param rate_law: reaction rate law(s) given as a string (?)
        :param enzymes: set of enzymes that catalyze the reaction
        """
        self.id = id
        self.name = name if name is not None else id
        self.reversible = reversible
        self.is_exchange = is_exchange
        self.stoichiometry = stoichiometry
        self.modifiers = modifiers
        self.parameters = parameters
        self.rate_law = rate_law
        self.enzymes = enzymes


class KineticModel:
    def __init__(self, model_id: str):
        """
        Constructor for representation of a system of metabolic reactions.

        All attributes apart from model_id are initialized as empty defaultdicts.

        Each of the dictionary will be of the form {'entity_id': entity_object}, where entity stands for metabolite,
        reaction, compartment, or condition, at the moment.
        """
        self.model_id = model_id
        self.metabolites = defaultdict()
        self.reactions = defaultdict()
        self.compartments = defaultdict()


class Measurement:
    def __init__(self,
                 target_id: str,
                 value: float,
                 uncertainty: float = None,
                 scale: str = None,
                 target_type: str = None):
        """
        Constructor for measurement object.

        :param target_id: id of the thing being measured
        :param value: value for the measurement
        :param uncertainty: uncertainty associated to the measurent
        :param scale: scale of the measurement, e.g. 'log10' or 'linear
        :param target_type: type of thing being measured, e.g. 'metabolite', 'reaction', 'enzyme'.
        """
        self.target_id = target_id
        self.value = value
        self.uncertainty = uncertainty
        self.scale = scale
        self.target_type = target_type


class Experiment:
    def __init__(self,
                 id: str,
                 met_meas: Dict[str, Measurement] = defaultdict(),
                 rxn_meas: Dict[str, Measurement] = defaultdict(),
                 enz_meas: Dict[str, Measurement] = defaultdict(),
                 metadata: str = None):

        """
        Constructor for condition object.

        :param id: condition id
        :param unbalanced_met_info:
        :param met_meas: metabolite measurements as a dictionary of the form {'met_id', measurement}
        :param rxn_meas: reaction measurements as a dictionary of the form {'rxn_id', measurement}
        :param enz_meas: enzyme measurements as a dictionary of the form {'enz_id', measurement}
        :param metadata: any info about the condition
        """
        self.id = id
        self.met_meas = met_meas
        self.rxn_meas = rxn_meas
        self.enz_meas = enz_meas
        self.metadata = metadata


class Prior:
    def __init__(self,
                 id: str,
                 target_id: str,
                 location: float,
                 scale: float,
                 target_type: str):
        """
        A prior distribuition.

        As currently implemented, the target must be a single parameter and the
        distribution must have a location and a scale.
        
        :param id: a string identifying the prior object
        :param target_id: a string identifying the thing that has a prior distribution.
        :param location: a number specifying the location of the distribution
        :param scale: a number specifying the scale of the distribution
        :param target_type: a string describing the target, e.g. 'kinetic_parameter' or 'unbalanced_metabolite'
        """
        self.id = id
        self.target_id = target_id
        self.location = location
        self.scale = scale
        self.target_type = target_type


class EnzymeKatInput:
    def __init__(self,
                 kinetic_model: KineticModel,
                 priors: Dict[str, Prior],
                 experiments: Dict[str, Experiment] = defaultdict()):
        """
        Everything that is needed to run EnzymeKat.
        
        :param kinetic_system: a KineticSystem object
        :param priors: a dictionary mapping prior ids to Prior objects
        :param experiments: a dictionary mapping experiment ids to Experiment objects
        """
        self.kinetic_model = kinetic_model
        self.priors = priors
        self.experiments = experiments
