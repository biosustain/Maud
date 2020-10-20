# Copyright (C) 2019 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Definitions of Maud-specific objects."""

from collections import defaultdict
from typing import Dict, List


class Compartment:
    """Constructor for compartment objects.

    :param id: compartment id, use a BiGG id if possible.
    :param name: compartment name.
    :param volume: compartment volume.
    """

    def __init__(self, id: str, name: str = None, volume: float = 1.0):
        self.id = id
        self.name = name if name is not None else id
        self.volume = volume


class Metabolite:
    """Constructor for metabolite objects.

    :param id: metabolite id, use a BiGG id if possible.
    :param name: metabolite name.
    """

    def __init__(self, id: str, name: str = None):
        self.id = id
        self.name = name if name is not None else id


class MetaboliteInCompartment:
    """A metabolite, in a compartment, or mic for short.

    :param id: this mic's id, usually <metabolite_id>_<compartment_id>.
    :param metabolite_id: id of this mic's metabolite
    :param compartment_id: id of this mic's compartment
    :param balanced: Does this mic have stable concentration at steady state?

    """

    def __init__(
        self,
        id: str,
        metabolite_id: str = None,
        compartment_id: str = None,
        balanced: bool = None,
    ):
        self.id = id
        self.metabolite_id = metabolite_id
        self.compartment_id = compartment_id
        self.balanced = balanced


class Modifier:
    """Constructor for modifier objects.

    :param mic_id: the id of the modifying metabolite-in-compartment
    :param enzyme_id: the id of the modified enzyme
    :param modifier_type: what is the modifier type, e.g.
    'allosteric_activator', 'allosteric_inhibitor', 'competitive inhibitor'
    """

    def __init__(self, mic_id: str, enzyme_id: str, modifier_type: str = None):
        allosteric_types = ["allosteric_inhibitor", "allosteric_activator"]
        self.mic_id = mic_id
        self.enzyme_id = enzyme_id
        self.modifier_type = modifier_type
        self.allosteric = modifier_type in allosteric_types


class Parameter:
    """Constructor for parameter object.

    :param id: parameter id
    :param enzyme_id: id of the enzyme associated with the parameter
    :param metabolite_id: id of the metabolite associated with the parameter if any
    """

    def __init__(
        self,
        id: str,
        enzyme_id: str,
        metabolite_id: str = None,
        is_thermodynamic: bool = False,
    ):
        self.id = id
        self.enzyme_id = enzyme_id
        self.metabolite_id = metabolite_id
        self.is_thermodynamic = is_thermodynamic


class Enzyme:
    """Constructor for the enzyme object.

    :param id: a string identifying the enzyme
    :param reaction_id: the id of the reaction the enzyme catalyses
    :param name: human-understandable name for the enzyme
    :param modifiers: modifiers, given as {'modifier_id': modifier_object}
    :param subunits: number of subunits in enzymes
    """

    def __init__(
        self,
        id: str,
        reaction_id: str,
        name: str,
        modifiers: Dict[str, List[Modifier]] = None,
        subunits: int = 1,
    ):
        if modifiers is None:
            modifiers = defaultdict()
        self.id = id
        self.name = name
        self.modifiers = modifiers
        self.subunits = subunits


class Reaction:
    """Constructor for the reaction object.

    :param id: reaction id, use a BiGG id if possible.
    :param name: reaction name.
    :param reversible: whether or not reaction is reversible.
    :param is_exchange: whether or not reaction is an exchange reaction.
    :param stoichiometry: reaction stoichiometry,
    e.g. for the reaction: 1.5 f6p <-> fdp we have {'f6p'; -1.5, 'fdp': 1}
    :param enzymes: Dictionary mapping enzyme ids to Enzyme objects
    :param water_stroichiometry: Reaction's stoichiometric coefficient for water
    """

    def __init__(
        self,
        id: str,
        name: str = None,
        reversible: bool = True,
        is_exchange: bool = None,
        stoichiometry: Dict[str, float] = None,
        enzymes: Dict[str, Enzyme] = None,
        water_stoichiometry: float = 0,
    ):
        if stoichiometry is None:
            stoichiometry = defaultdict()
        if enzymes is None:
            enzymes = defaultdict()
        self.id = id
        self.name = name if name is not None else id
        self.reversible = reversible
        self.is_exchange = is_exchange
        self.stoichiometry = stoichiometry
        self.enzymes = enzymes
        self.water_stoichiometry = water_stoichiometry

class Drain:
    """Constructor for the reaction object.

    :param id: drain id, use a BiGG id if possible.
    :param name: drain name.
    :param stoichiometry: reaction stoichiometry,
    """

    def __init__(
        self,
        id: str,
        name: str = None,
        stoichiometry: Dict[str, float] = None,
    ):
        if stoichiometry is None:
            stoichiometry = defaultdict()
        if enzymes is None:
            enzymes = defaultdict()
        self.id = id
        self.name = name if name is not None else id
        self.stoichiometry = stoichiometry

class KineticModel:
    """Constructor for representation of a system of metabolic reactions.

    :param model_id: id of the kinetic model
    :param metabolites: dictionary mapping strings to metabolite objects
    :param reactions: dictionary mapping strings to reaction objects
    :param compartments: dictionary mapping strings to compartment objects
    :param mic: dictionary mapping strings to MetaboliteInCompartment objects
    """

    def __init__(
        self,
        model_id: str,
        metabolites: Dict[str, Metabolite],
        reactions: Dict[str, Reaction],
        drains: Dict[str, Drain],
        compartments: Dict[str, Compartment],
        mics: Dict[str, MetaboliteInCompartment],
    ):
        self.model_id = model_id
        self.metabolites = metabolites
        self.reactions = reactions
        self.drains = drains
        self.compartments = compartments
        self.mics = mics


class Measurement:
    """Constructor for measurement object.

    :param target_id: id of the thing being measured
    :param value: value for the measurement
    :param uncertainty: uncertainty associated to the measurent
    :param scale: scale of the measurement, e.g. 'log10' or 'linear
    :param target_type: type of thing being measured, e.g. 'metabolite', 'reaction',
    'enzyme', 'drain'.
    """

    def __init__(
        self,
        target_id: str,
        value: float,
        uncertainty: float = None,
        scale: str = None,
        target_type: str = None,
    ):
        self.target_id = target_id
        self.value = value
        self.uncertainty = uncertainty
        self.scale = scale
        self.target_type = target_type


class Experiment:
    """Constructor for condition object.

    :param id: condition id
    :param unbalanced_met_info:
    :param measurements: dictionary mapping keys 'enzyme', 'metabolite' and 'reaction'
    to dictionaries with the form {target id: measurement}
    :param metadata: any info about the condition
    """

    def __init__(
        self,
        id: str,
        measurements: Dict[str, Dict[str, Measurement]] = None,
        metadata: str = None,
        knockouts: List[str] = None,
    ):
        if measurements is None:
            measurements = defaultdict()
        self.id = id
        self.measurements = measurements
        self.metadata = metadata
        self.knockouts = knockouts


class Prior:
    """A prior distribuition.

    As currently implemented, the target must be a single parameter and the
    distribution must have a location and a scale.

    :param id: a string identifying the prior object
    :param target_id: a string identifying the thing that has a prior distribution
    :param location: a number specifying the location of the distribution
    :param scale: a number specifying the scale of the distribution
    :param experiment_id: id of the relevant experiment (for enzymes or unbalanced
    metabolites)
    :param mic_id: id of relevant metabolite-in-compartment
    :param metabolite_id: id of relevant metabolite
    :param enzyme_id: id of relevant enzyme
    :param drain_id: id of relevant drain
    """

    def __init__(
        self,
        id: str,
        location: float,
        scale: float,
        experiment_id: str = None,
        mic_id: str = None,
        metabolite_id: str = None,
        enzyme_id: str = None,
        drain_id: str = None,
    ):
        self.id = id
        self.location = location
        self.scale = scale
        self.experiment_id = experiment_id
        self.mic_id = mic_id
        self.metabolite_id = metabolite_id
        self.enzyme_id = enzyme_id


class MaudInput:
    """Everything that is needed to run Maud.

    :param kinetic_system: a KineticSystem object
    :param priors: a dictionary mapping prior types to lists of Prior objects
    :param stan codes: a dictionary with keys 'metabolite', 'reaction', 'mic',
    'experiment', 'balanced_mic', 'unbalanced_mic' and 'kinetic_parameter', whose
    values are dictionaries mapping ids of the relevant objects to integer codes
    :param experiments: a dictionary mapping experiment ids to Experiment objects
    """

    def __init__(
        self,
        kinetic_model: KineticModel,
        priors: Dict[str, List[Prior]],
        stan_codes: Dict[str, Dict[str, int]],
        experiments: Dict[str, Experiment] = None,
    ):
        if experiments is None:
            experiments = defaultdict()
        self.kinetic_model = kinetic_model
        self.priors = priors
        self.stan_codes = stan_codes
        self.experiments = experiments
