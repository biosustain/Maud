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
from dataclasses import dataclass
from typing import Dict, List

import numpy as np
from cmdstanpy import CmdStanMCMC

from maud.utils import (
    get_lognormal_parameters_from_quantiles,
    get_normal_parameters_from_quantiles,
)


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
    :param water_stroichiometry: Reaction's stoichiometric coefficient for water
    """

    def __init__(
        self,
        id: str,
        reaction_id: str,
        name: str,
        modifiers: Dict[str, List[Modifier]] = None,
        subunits: int = 1,
        water_stoichiometry: float = 0,
    ):
        if modifiers is None:
            modifiers = defaultdict()
        self.id = id
        self.name = name
        self.modifiers = modifiers
        self.subunits = subunits
        self.water_stoichiometry = water_stoichiometry


class Reaction:
    """Constructor for the reaction object.

    :param id: reaction id, use a BiGG id if possible.
    :param name: reaction name.
    :param reversible: whether or not reaction is reversible.
    :param is_exchange: whether or not reaction is an exchange reaction.
    :param stoichiometry: reaction stoichiometry,
    e.g. for the reaction: 1.5 f6p <-> fdp we have {'f6p'; -1.5, 'fdp': 1}
    :param enzymes: Dictionary mapping enzyme ids to Enzyme objects
    """

    def __init__(
        self,
        id: str,
        name: str,
        stoichiometry: Dict[str, float],
        enzymes: List[Enzyme],
    ):
        if stoichiometry is None:
            stoichiometry = defaultdict()
        if enzymes is None:
            enzymes = defaultdict()
        self.id = id
        self.name = name if name is not None else id
        self.stoichiometry = stoichiometry
        self.enzymes = enzymes


class Drain:
    """Constructor for the reaction object.

    :param id: drain id, use a BiGG id if possible.
    :param name: drain name.
    :param stoichiometry: reaction stoichiometry,
    """

    def __init__(
        self, id: str, name: str = None, stoichiometry: Dict[str, float] = None
    ):
        if stoichiometry is None:
            stoichiometry = defaultdict()
        self.id = id
        self.name = name if name is not None else id
        self.stoichiometry = stoichiometry


class Phosphorylation:
    """Constructor for the phosphorylation object.

    :param id: phosphorylation id. use BIGG id if possible.
    :param name: name of phosphorylation reaction.
    :param activating: if the interaction activates the
    target enzyme.
    :param inhibiting: if the interaction inhibits the
    target enzyme.
    :enzyme_id: the target enzyme of the interaction
    """

    def __init__(
        self,
        id: str,
        name: str = None,
        activating: bool = None,
        inhibiting: bool = None,
        enzyme_id: str = None,
    ):
        self.id = id
        self.name = name
        self.activating = activating
        self.inhibiting = inhibiting
        self.enzyme_id = enzyme_id


class KineticModel:
    """Constructor for representation of a system of metabolic reactions.

    :param model_id: id of the kinetic model
    :param metabolites: list of metabolite objects
    :param reactions: list of reaction objects
    :param drains: list of drain objects
    :param compartments: list of compartment objects
    :param mic: list of MetaboliteInCompartment objects
    """

    def __init__(
        self,
        model_id: str,
        metabolites: List[Metabolite],
        reactions: List[Reaction],
        compartments: List[Compartment],
        mics: List[MetaboliteInCompartment],
        drains: List[Drain] = None,
        phosphorylation: List[Phosphorylation] = None,
    ):
        self.model_id = model_id
        self.metabolites = metabolites
        self.reactions = reactions
        self.drains = drains if drains is not None else []
        self.compartments = compartments
        self.mics = mics
        self.phosphorylation = phosphorylation if phosphorylation is not None else []


class Measurement:
    """Constructor for measurement object.

    :param target_id: id of the thing being measured
    :param experiment_id: id of the experiment where the measurement was done
    :param value: measured value
    :param error: error associated to the measurent
    :param target_type: type of thing being measured, e.g. 'metabolite', 'reaction',
    'enzyme'.
    """

    def __init__(
        self,
        target_id: str,
        experiment_id: str,
        target_type: str,
        value: float,
        error: float,
    ):
        self.target_id = target_id
        self.experiment_id = experiment_id
        self.target_type = target_type
        self.value = value
        self.error = error


class Knockout:
    """Constructor for knockout object.

    :param experiment_id: id of the experiment where the thing was knocked out
    :target_id: id of the thing that was knocked out
    :knockout_type: either "enz" or "phos"
    """

    def __init__(
        self,
        experiment_id: str,
        target_id: str,
        knockout_type: str,
    ):
        self.experiment_id = experiment_id
        self.target_id = target_id
        self.knockout_type = knockout_type


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
    :param phos_enz_id: id of enzyme involved with phosphorylation reaction
    :param drain_id: id of relevant drain
    """

    def __init__(
        self,
        id: str,
        is_non_negative: bool,
        experiment_id: str = None,
        mic_id: str = None,
        metabolite_id: str = None,
        enzyme_id: str = None,
        phos_enz_id: str = None,
        drain_id: str = None,
        pct1: float = None,
        pct99: float = None,
        location: float = None,
        scale: float = None,
    ):
        self.id = id
        self.is_non_negative = is_non_negative
        self.experiment_id = experiment_id
        self.mic_id = mic_id
        self.metabolite_id = metabolite_id
        self.enzyme_id = enzyme_id
        self.phos_enz_id = phos_enz_id
        self.drain_id = drain_id

        self.pct1 = pct1
        self.pct99 = pct99
        if pct1 is not None and pct99 is not None:
            if is_non_negative:
                mu, self.scale = get_lognormal_parameters_from_quantiles(
                    pct1, 0.01, pct99, 0.99
                )
                self.location = np.exp(mu)
            else:
                (
                    self.location,
                    self.scale,
                ) = get_normal_parameters_from_quantiles(pct1, 0.01, pct99, 0.99)
        else:
            self.location = location
            self.scale = scale


class PriorSet:
    """Object containing all priors for a MaudInput."""

    def __init__(
        self,
        kcat_priors: List[Prior],
        phos_kcat_priors: List[Prior],
        km_priors: List[Prior],
        formation_energy_priors: List[Prior],
        unbalanced_metabolite_priors: List[Prior],
        inhibition_constant_priors: List[Prior],
        tense_dissociation_constant_priors: List[Prior],
        relaxed_dissociation_constant_priors: List[Prior],
        transfer_constant_priors: List[Prior],
        drain_priors: List[Prior],
        enzyme_concentration_priors: List[Prior],
        phos_enz_concentration_priors: List[Prior],
    ):
        self.kcat_priors = kcat_priors
        self.phos_kcat_priors = phos_kcat_priors
        self.km_priors = km_priors
        self.formation_energy_priors = formation_energy_priors
        self.unbalanced_metabolite_priors = unbalanced_metabolite_priors
        self.inhibition_constant_priors = inhibition_constant_priors
        self.tense_dissociation_constant_priors = tense_dissociation_constant_priors
        self.relaxed_dissociation_constant_priors = relaxed_dissociation_constant_priors
        self.transfer_constant_priors = transfer_constant_priors
        self.drain_priors = drain_priors
        self.enzyme_concentration_priors = enzyme_concentration_priors
        self.phos_enz_concentration_priors = phos_enz_concentration_priors


@dataclass
class StanCoordSet:
    """Object containing human-readable indexes for Maud's parameters.

    These are "coordinates" in the sense of xarray

    """

    metabolites: List[str]
    mics: List[str]
    balanced_mics: List[str]
    unbalanced_mics: List[str]
    kms: List[str]
    reactions: List[str]
    experiments: List[str]
    enzymes: List[str]
    drains: List[str]
    phos_enzs: List[str]
    yconc_exps: List[str]
    yconc_mics: List[str]
    yflux_exps: List[str]
    yflux_rxns: List[str]
    yenz_exps: List[str]
    yenz_enzs: List[str]
    ci_mics: List[str]
    ai_mics: List[str]
    aa_mics: List[str]


class MaudConfig:
    """User's configuration for a Maud input.

    :param name: name for the input. Used to name the output directory
    :param kinetic_model_file: path to a valid kientic model file.
    :param priors_file: path to a valid priors file.
    :param experiments_file: path to a valid experiments file.
    :param likelihood: Whether or not to take measurements into account.
    :param ode_config: Configuration for Stan's ode solver.
    :param cmdstanpy_config: Arguments to cmdstanpy.CmdStanModel.sample.
    """

    def __init__(
        self,
        name: str,
        kinetic_model_file: str,
        priors_file: str,
        experiments_file: str,
        likelihood: bool,
        ode_config: dict,
        cmdstanpy_config: dict,
    ):
        self.name = name
        self.kinetic_model_file = kinetic_model_file
        self.priors_file = priors_file
        self.experiments_file = experiments_file
        self.likelihood = likelihood
        self.ode_config = ode_config
        self.cmdstanpy_config = cmdstanpy_config


class MaudInput:
    """Everything that is needed to run Maud.

    :param kinetic_system: a KineticSystem object
    :param priors: a dictionary mapping prior types to lists of Prior objects
    :param stan codes: a dictionary with keys 'metabolite', 'reaction', 'mic',
    'experiment', 'balanced_mic', 'unbalanced_mic' and 'kinetic_parameter', whose
    values are dictionaries mapping ids of the relevant objects to integer codes
    :param measurement_set: a list of Measurement objects
    :param knockout_set: a list of Knockout objects
    """

    def __init__(
        self,
        config: MaudConfig,
        kinetic_model: KineticModel,
        priors: PriorSet,
        stan_coords: StanCoordSet,
        measurements: List[Measurement],
        knockouts: List[Knockout],
    ):
        self.config = config
        self.kinetic_model = kinetic_model
        self.priors = priors
        self.stan_coords = stan_coords
        self.measurements = measurements
        self.knockouts = knockouts


class SimulationStudyOutput:
    """Expected output of a simulation study.

    :param input_data_sim: dictionary used to create simulation
    :param input_data_sample: dictionary used to create samples
    :param true_values: dictionary mapping param names to true values
    :param sim: CmdStanMCMC that generated simulated measurements
    :param mi: Maud input used for sampling
    :param samples: CmdStanMCMC output of the simulation study
    """

    def __init__(
        self,
        input_data_sim: dict,
        input_data_sample: dict,
        true_values: dict,
        sim: CmdStanMCMC,
        mi: MaudInput,
        samples: CmdStanMCMC,
    ):
        self.input_data_sim = input_data_sim
        self.input_data_sample = input_data_sample
        self.true_values = true_values
        self.sim = sim
        self.mi = mi
        self.samples = samples
