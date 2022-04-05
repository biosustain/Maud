# Copyright (C) 2022 Novo Nordisk Foundation Center for Biosustainability,
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

"""Maud's definition of a kinetic model."""

from dataclasses import field
from enum import Enum
from itertools import chain
from typing import Dict, List

from pydantic import BaseModel, validator
from pydantic.dataclasses import dataclass


ID_SEPARATOR = "_"


class ReactionMechanism(str, Enum):
    RMM = "reversible_michaelis_menten"
    IMM = "irreversible_michaelis_menten"
    DRAIN = "drain"


class ModificationType(str, Enum):
    ACTIVATION = "activation"
    INHIBITION = "inhibition"


@dataclass
class Metabolite(BaseModel):
    id: str
    name: str
    inchi_key: str

    @validator("id")
    def id_must_not_contain_seps(cls, v):
        assert ID_SEPARATOR not in v
        return v


@dataclass
class Enzyme(BaseModel):
    id: str
    name: str
    subunits: int

    @validator("id")
    def id_must_not_contain_seps(cls, v):
        assert ID_SEPARATOR not in v
        return v

    @validator("subunits")
    def subunits_must_be_positive(cls, v):
        assert v > 0
        return v


@dataclass
class Compartment(BaseModel):
    id: str
    name: str
    volume: float

    @validator("id")
    def id_must_not_contain_seps(cls, v):
        assert ID_SEPARATOR not in v
        return v


@dataclass
class Reaction:
    id: str
    name: str
    mechanism: ReactionMechanism
    stoichiometry: Dict[str, float]
    water_stoichiometry: float

    @validator("id")
    def id_must_not_contain_seps(cls, v):
        assert ID_SEPARATOR not in v
        return v

    @validator("stoichiometry")
    def stoichiometry_must_be_non_zero(cls, v):
        assert v != 0
        return v


@dataclass
class MetaboliteInCompartment:
    metabolite_id: str
    compartment_id: str
    balanced: bool
    id: str = field(init=False)

    def __post_init__(self):
        self.id = self.metabolite_id + ID_SEPARATOR + self.compartment_id


@dataclass
class EnzymeReaction:
    enzyme_id: str
    reaction_id: str
    id: str = field(init=False)

    def __post_init__(self):
        self.id = self.enzyme_id + ID_SEPARATOR + self.reaction_id


@dataclass
class Allostery:
    enzyme_id: str
    mic_id: str
    modification_type: ModificationType
    id: str = field(init=False)

    def __post_init__(self):
        self.id = (
            self.enzyme_id
            + ID_SEPARATOR
            + self.mic_id
            + ID_SEPARATOR
            + self.modification_type
        )


@dataclass
class CompetitiveInhibition:
    er_id: str
    mic_id: str
    modification_type: ModificationType
    id: str = field(init=False)

    def __post_init__(self):
        self.id = (
            self.er_id
            + ID_SEPARATOR
            + self.mic_id
            + ID_SEPARATOR
            + self.modification_type
        )


@dataclass
class Phosphorylation:
    enzyme_id: str
    modification_type: ModificationType
    id: str = field(init=False)

    def __post_init__(self):
        self.id = self.enzyme_id + ID_SEPARATOR + self.modification_type


class KineticModel(BaseModel):
    """Representation of a system of metabolic network."""

    name: str
    metabolites: List[Metabolite]
    enzymes: List[Enzyme]
    compartments: List[Compartment]
    reactions: List[Reaction]
    mics: List[MetaboliteInCompartment]
    ers: List[EnzymeReaction]
    allosteries: List[Allostery]
    competitive_inhibitions: List[CompetitiveInhibition]
    phosphorylations: List[Phosphorylation]

    @validator("metabolites")
    def metabolite_ids_must_be_unique(cls, v):
        met_ids = [m.id for m in v]
        assert len(met_ids) == len(set(met_ids))
        return v

    @validator("enzymes")
    def enzyme_ids_must_be_unique(cls, v):
        met_ids = [m.id for m in v]
        assert len(met_ids) == len(set(met_ids))
        return v

    @validator("compartments")
    def compartment_ids_must_be_unique(cls, v):
        met_ids = [m.id for m in v]
        assert len(met_ids) == len(set(met_ids))
        return v

    @validator("reactions")
    def reaction_ids_must_be_unique(cls, v):
        rxn_ids = [r.id for r in v]
        assert len(rxn_ids) == len(set(rxn_ids))
        return v

    @validator("reactions")
    def stoic_keys_must_be_mic_ids(cls, v, values):
        for k in chain.from_iterable(r.stoichiometry.keys() for r in v):
            assert k in values["mics"]
        return v

    @validator("mics")
    def mic_references_must_exist(cls, v, values):
        metabolite_ids = [m.id for m in values["metabolites"]]
        compartment_ids = [m.id for m in values["compartments"]]
        for mic in v:
            assert mic.metabolite_id in metabolite_ids
            assert mic.compartment_id in compartment_ids
        return v

    @validator("ers")
    def er_references_must_exist(cls, v, values):
        enzyme_ids = [m.id for m in values["enzymes"]]
        reaction_ids = [m.id for m in values["reactions"]]
        for er in v:
            assert er.enzyme_id in enzyme_ids
            assert er.reaction_id in reaction_ids
        return v

    @validator("allosteries")
    def allostery_references_must_exist(cls, v, values):
        enzyme_ids = [m.id for m in values["enzymes"]]
        mic_ids = [m.id for m in values["mics"]]
        for allostery in v:
            assert allostery.enzyme_id in enzyme_ids
            assert allostery.mic_id in mic_ids
        return v

    @validator("competitive_inhibitions")
    def ci_references_must_exist(cls, v, values):
        er_ids = [m.id for m in values["ers"]]
        mic_ids = [m.id for m in values["mics"]]
        for ci in v:
            assert ci.enzyme_id in er_ids
            assert ci.mic_id in mic_ids
        return v

    @validator("phosphorylations")
    def phosphorylation_references_must_exist(cls, v, values):
        enzyme_ids = [m.id for m in values["enzymes"]]
        for ci in v:
            assert ci.enzyme_id in enzyme_ids
        return v
