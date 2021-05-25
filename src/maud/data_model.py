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
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from cmdstanpy import CmdStanMCMC


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
        self.reaction_id = reaction_id
        self.name = name
        self.modifiers = modifiers
        self.subunits = subunits
        self.water_stoichiometry = water_stoichiometry
        self.allosteric = (
            len(self.modifiers["allosteric_activator"]) > 0
            or len(self.modifiers["allosteric_inhibitor"]) > 0
        )


class Reaction:
    """Constructor for the reaction object.

    :param id: reaction id, use a BiGG id if possible.
    :param name: reaction name.
    :param reaction_type: either "reversible_modular_rate_law" or "drain".
    :param stoichiometry: reaction stoichiometry,
    e.g. for the reaction: 1.5 f6p <-> fdp we have {'f6p'; -1.5, 'fdp': 1}
    :param enzymes: Dictionary mapping enzyme ids to Enzyme objects
    """

    def __init__(
        self,
        id: str,
        name: str,
        reaction_type: str,
        stoichiometry: Dict[str, float],
        enzymes: List[Enzyme],
    ):
        if stoichiometry is None:
            stoichiometry = defaultdict()
        if enzymes is None:
            enzymes = defaultdict()
        self.id = id
        self.name = name if name is not None else id
        self.reaction_type = reaction_type
        self.stoichiometry = stoichiometry
        self.enzymes = enzymes


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
        phosphorylation: List[Phosphorylation] = None,
    ):
        self.model_id = model_id
        self.metabolites = metabolites
        self.reactions = reactions
        self.compartments = compartments
        self.mics = mics
        self.phosphorylation = phosphorylation if phosphorylation is not None else []


@dataclass
class MeasurementSet:
    """A container for a complete set of measurements, including knockouts."""

    yconc: pd.DataFrame
    yflux: pd.DataFrame
    yenz: pd.DataFrame
    enz_knockouts: pd.DataFrame
    phos_knockouts: pd.DataFrame


@dataclass
class IndPrior1d:
    """Independent location/scale prior for a 1-dimentional parameter."""

    parameter_name: str
    location: pd.Series
    scale: pd.Series

    def __post_init__(self):
        """Check that location and scale indexes match."""
        if not self.location.index.equals(self.scale.index):
            raise ValueError("Location index doesn't match scale index.")


@dataclass
class IndPrior2d:
    """Independent location/scale prior for a 2-dimensional parameter."""

    parameter_name: str
    location: pd.DataFrame
    scale: pd.DataFrame

    def __post_init__(self):
        """Check that location and scale indexes and columns match."""
        if not self.location.index.equals(self.scale.index):
            raise ValueError("Location index doesn't match scale index.")
        if not self.location.columns.equals(self.scale.columns):
            raise ValueError("Location columns don't match scale columns.")


@dataclass
class PriorSet:
    """Object containing all priors for a MaudInput."""

    priors_kcat: IndPrior1d
    priors_kcat_phos: IndPrior1d
    priors_km: IndPrior1d
    priors_dgf: IndPrior1d
    priors_ki: IndPrior1d
    priors_diss_t: IndPrior1d
    priors_diss_r: IndPrior1d
    priors_transfer_constant: IndPrior1d
    priors_conc_unbalanced: IndPrior2d
    priors_drain: IndPrior2d
    priors_conc_enzyme: IndPrior2d
    priors_conc_phos: IndPrior2d


@dataclass
class StanCoordSet:
    """Object containing human-readable indexes for Maud's parameters.

    These are "coordinates" in the sense of xarray

    """

    metabolites: List[str]
    mics: List[str]
    balanced_mics: List[str]
    unbalanced_mics: List[str]
    km_enzs: List[str]
    km_mics: List[str]
    reactions: List[str]
    experiments: List[str]
    enzymes: List[str]
    edges: List[str]
    allosteric_enzymes: List[str]
    drains: List[str]
    phos_enzs: List[str]
    yconc_exps: List[str]
    yconc_mics: List[str]
    yflux_exps: List[str]
    yflux_rxns: List[str]
    yenz_exps: List[str]
    yenz_enzs: List[str]
    ci_enzs: List[str]
    ci_mics: List[str]
    ai_enzs: List[str]
    ai_mics: List[str]
    aa_enzs: List[str]
    aa_mics: List[str]
    enz_ko_exps: List[str]
    enz_ko_enzs: List[str]
    phos_ko_exps: List[str]
    phos_ko_enzs: List[str]


class MaudConfig:
    """User's configuration for a Maud input.

    :param name: name for the input. Used to name the output directory
    :param kinetic_model_file: path to a valid kientic model file.
    :param priors_file: path to a valid priors file.
    :param experiments_file: path to a valid experiments file.
    :param likelihood: Whether or not to take measurements into account.
    :param ode_config: Configuration for Stan's ode solver.
    :param cmdstanpy_config: Arguments to cmdstanpy.CmdStanModel.sample.
    :param user_inits_file: path to a csv file of initial values.
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
        user_inits_file: Optional[str],
    ):
        self.name = name
        self.kinetic_model_file = kinetic_model_file
        self.priors_file = priors_file
        self.experiments_file = experiments_file
        self.likelihood = likelihood
        self.ode_config = ode_config
        self.cmdstanpy_config = cmdstanpy_config
        self.user_inits_file = user_inits_file


class MaudInput:
    """Everything that is needed to run Maud.

    :param kinetic_system: a KineticSystem object
    :param priors: a dictionary mapping prior types to lists of Prior objects
    :param stan_coords: a StanCoordSet object
    :param measurement_set: a list of Measurement objects
    :param inits: a dictionary of initial parameter values
    """

    def __init__(
        self,
        config: MaudConfig,
        kinetic_model: KineticModel,
        priors: PriorSet,
        stan_coords: StanCoordSet,
        measurements: MeasurementSet,
        inits: Dict[str, np.array],
    ):
        self.config = config
        self.kinetic_model = kinetic_model
        self.priors = priors
        self.stan_coords = stan_coords
        self.measurements = measurements
        self.inits = inits


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
