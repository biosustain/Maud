from collections import defaultdict
import toml
from typing import Dict
​
​
class Compartment:
    def __init__(self, id : str, name: str = None, volume: float = 1.0):
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
​
​
class Metabolite:
    def __init__(self, id: str, name: str = None, balanced: bool = None, compartment: Compartment = None):
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
​
​
class Reaction:
    def __init__(self, id: str, name: str = None, reversible: bool = True, is_exchange: bool = None,
                 stoichiometry: dict = None, modifiers: dict = None, parameters: dict = None, rate_law: str = None,
                 enzymes: set = None):
        """
        Constructor for the reaction object.
​
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


class Measurement:
    def __init__(self, value: float, uncertainty: float = None, scale: str = None, measurement_type: str = None):
        """
        Constructor for measurement object, assuming option 2 for the condition class.
​
        :param value: value for the measurement
        :param uncertainty: uncertainty associated to the measurent
        :param scale: scale of the measurement, e.g. 'log10' or 'linear
        :param measurement_type: type of the measurement, e.g. 'met concentration', 'rxn flux', 'enzyme concentration', etc.
        """
        self.value = value
        self.uncertainty = uncertainty
        self.scale = scale
        self.measurement_type = measurement_type
​
​
class Experiment:
    def __init__(self,
                 id: str,
                 met_meas: Dict[str, Measurement] = None,
                 rxn_meas: Dict[str, Measurement] = None,
                 enz_meas: Dict[str, Measurement] = None,
                 metadata: str = None):
​
        """
        Constructor for condition object.
​
        :param id: condition id
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
​
​
class EnzymeKatInput:
    def __init__(self, model_id):
        """
        Constructor for the model object,
        All attributes apart from model_id are initialized as empty defaultdics.
        Each of the dictionary will be of the form {'entity_id': entity_object}, where entity stands for metabolite,
        reaction, compartment, or condition, at the moment.
        """
        self.model_id = model_id
        self.metabolites = defaultdict()
        self.reactions = defaultdict()
        self.compartments = defaultdict()
        self.experiments = defaultdict()

    
def load_toml(filepath: str, id: str ='eki') -> EnzymeKatInput:
    """
    Load an EnzymeKatInput object from a suitable toml file
    
    :param filepath: path to a toml file
    :param id: id for the output object

    """
    eki = EnzymeKatInput(id)
    parsed_toml = toml.load(filepath)
    for c in parsed_toml['compartments']:
        cmp = Compartment(
            id=c['id'],
            name=c['name'],
            volume=c['volume']
        )
        eki.compartments.update({cmp.id: cmp})
    for m in parsed_toml['metabolites']:
        met = Metabolite(
            id=m['id'],
            name=m['name'],
            balanced=m['balanced'],
            compartment=m['compartment']
        )
        eki.metabolites.update({met.id: met})
    for r in parsed_toml['reactions']:
        rxn = Reaction(
            id=r['id'],
            name=r['name'] if r['name'] else None,
            reversible=r['reversible'] if 'reversible' in r.keys() else None,
            is_exchange=r['is_exchange'] if 'is_exchange' in r.keys() else None,
            stoichiometry=r['stoichiometry'],
            modifiers=r['modifiers'] if 'modifiers' in r.keys() else None,
            parameters=r['parameters'] if 'parameters' in r.keys() else None,
            rate_law=r['rate_law'] if 'rate_law' in r.keys() else None,
            enzymes=r['enzymes'] if 'enzymes' in r.keys() else None 
        )
        eki.reactions.update({rxn.id: rxn})
    for e in parsed_toml['experiments']:
        measurements = {}
        for measurement_type in ['metabolite', 'enzyme', 'reaction']:
            type_msmts = {}
            for m in e['measurements'][measurement_type]:
                msmt = Measurement(
                    value=m['value'],
                    uncertainty=m['uncertainty'],
                    scale='ln',
                    measurement_type=measurement_type
                )
                type_msmts.update({m['target']: msmt})
            measurements[measurement_type] = type_msmts
        exp = Experiment(
            id=e['id'],
            met_meas=measurements['metabolite'],
            rxn_meas=measurements['reaction'],
            enz_meas=measurements['enzyme'],
            metadata=e['metadata']
        )
        eki.experiments.update({exp.id: exp})
    return eki

