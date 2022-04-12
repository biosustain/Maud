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
from typing import Dict, List, Optional, Union

import pandas as pd
from pydantic import root_validator, validator
# hack to use normal dataclass when type checking
# see https://github.com/microsoft/pyright/issues/1510
from pydantic.dataclasses import dataclass

from maud.data_model.hardcoding import ID_SEPARATOR


class ReactionMechanism(str, Enum):
    REVERSIBLE_MICHAELIS_MENTEN = 1
    IRREVERSIBLE_MICHAELIS_MENTEN = 2
    DRAIN = 3


class ModificationType(str, Enum):
    ACTIVATION = 1
    INHIBITION = 2


class KMConfig:
    arbitrary_types_allowed = True


@dataclass
class Metabolite:
    id: str
    name: Optional[str]
    inchi_key: Optional[str]

    @validator("id")
    def id_must_not_contain_seps(cls, v):
        assert ID_SEPARATOR not in v
        return v


@dataclass
class Enzyme:
    id: str
    name: Optional[str]
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
class Compartment:
    id: str
    name: Optional[str]
    volume: float

    @validator("id")
    def id_must_not_contain_seps(cls, v):
        assert ID_SEPARATOR not in v
        return v


@dataclass
class Reaction:
    id: str
    name: Optional[str]
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
        self.id = ID_SEPARATOR.join([self.metabolite_id, self.compartment_id])


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
    metabolite_id: str
    compartment_id: str
    modification_type: ModificationType
    id: str = field(init=False)

    def __post_init__(self):
        self.id = ID_SEPARATOR.join(
            [
                self.enzyme_id,
                self.metabolite_id,
                self.compartment_id,
            ]
        )
        self.mic_id = ID_SEPARATOR.join(
            [self.metabolite_id, self.compartment_id]
        )


@dataclass
class CompetitiveInhibition:
    enzyme_id: str
    reaction_id: str
    metabolite_id: str
    compartment_id: str
    id: str = field(init=False)

    def __post_init__(self):
        self.id = ID_SEPARATOR.join(
            [
                self.enzyme_id,
                self.reaction_id,
                self.metabolite_id,
                self.compartment_id,
            ]
        )
        self.er_id = ID_SEPARATOR.join([self.enzyme_id, self.reaction_id])
        self.mic_id = ID_SEPARATOR.join(
            [self.metabolite_id, self.compartment_id]
        )


@dataclass
class Phosphorylation:
    enzyme_id: str
    modification_type: ModificationType
    id: str = field(init=False)

    def __post_init__(self):
        self.id = self.enzyme_id + ID_SEPARATOR + self.modification_type


@dataclass(config=KMConfig)
class KineticModel:
    """Representation of a system of metabolic network."""

    name: str
    metabolites: List[Metabolite]
    enzymes: List[Enzyme]
    compartments: List[Compartment]
    mics: List[MetaboliteInCompartment]
    reactions: List[Reaction]
    ers: List[EnzymeReaction]
    allosteries: Optional[List[Allostery]]
    competitive_inhibitions: Optional[List[CompetitiveInhibition]]
    phosphorylations: Optional[List[Phosphorylation]]
    drains: List[Reaction] = field(init=False)
    edges: List[Union[Reaction, EnzymeReaction]] = field(init=False)
    stoichiometric_matrix: pd.DataFrame = field(init=False)

    def __post_init__(self):
        self.drains = [
            r for r in self.reactions if r.mechanism == ReactionMechanism.DRAIN
        ]
        self.edges = self.drains + self.ers
        self.stoichiometric_matrix = (
            get_stoichiometric_matrix(self.edges, self.mics, self.reactions)
        )

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

    @root_validator
    def stoic_keys_must_be_mic_ids(cls, values):
        mic_ids = [mic.id for mic in values["mics"]]
        for r in values["reactions"]:
            for stoich_mic_id in r.stoichiometry.keys():
                assert (
                    stoich_mic_id in mic_ids
                ), f"{r.id} has stoichiometry for bad mic_id {stoich_mic_id}"
        return values

    @root_validator
    def mic_references_must_exist(cls, values):
        metabolite_ids = [m.id for m in values["metabolites"]]
        compartment_ids = [c.id for c in values["compartments"]]
        for mic in values["mics"]:
            assert (
                mic.metabolite_id in metabolite_ids
            ), f"{mic.id} has bad metabolite_id."
            assert (
                mic.compartment_id in compartment_ids
            ), f"{mic.id} has bad compartment_id."
        return values

    @root_validator
    def er_references_must_exist(cls, values):
        enzyme_ids = [e.id for e in values["enzymes"]]
        reaction_ids = [r.id for r in values["reactions"]]
        for er in values["ers"]:
            assert er.enzyme_id in enzyme_ids, f"{er.id} has bad enzyme_id"
            assert (
                er.reaction_id in reaction_ids
            ), f"{er.id} has bad reaction_id"
        return values

    @root_validator
    def allostery_references_must_exist(cls, values):
        if values["allosteries"] is None:
            return values
        enzyme_ids = [e.id for e in values["enzymes"]]
        mic_ids = [mic.id for mic in values["mics"]]
        for allostery in values["allosteries"]:
            assert (
                allostery.enzyme_id in enzyme_ids
            ), f"{allostery.id} has bad enzyme_id"
            assert (
                allostery.mic_id in mic_ids
            ), f"{allostery.id} has bad mic_id"
        return values

    @root_validator
    def ci_references_must_exist(cls, values):
        if values["competitive_inhibitions"] is None:
            return values
        er_ids = [er.id for er in values["ers"]]
        mic_ids = [mic.id for mic in values["mics"]]
        for ci in values["competitive_inhibitions"]:
            assert ci.er_id in er_ids, f"{ci.id} has bad er_id"
            assert ci.mic_id in mic_ids, f"{ci.mic_id} has bad mic_id"
        return values

    @root_validator
    def phosphorylation_references_must_exist(cls, values):
        if values["phosphorylations"] is None:
            return values
        enzyme_ids = [e.id for e in values["enzymes"]]
        for ci in values["phosphorylations"]:
            assert ci.enzyme_id in enzyme_ids, f"{ci.id} has bad enzyme_id"
        return values


def get_stoichiometric_matrix(
    edges: List[Union[Reaction, EnzymeReaction]],
    mics: List[MetaboliteInCompartment],
    rxns: List[Reaction]
) -> pd.DataFrame:
    edge_ids = [e.id for e in edges]
    mic_ids = [mic.id for mic in mics]
    S = pd.DataFrame(0, index=mic_ids, columns=edge_ids)
    for e in edges:
        if isinstance(e, Reaction):
            rxn = e
        else:
            rxn = next(r for r in rxns if r.id == e.reaction_id)
        for mic_id, stoic in rxn.stoichiometry.items():
            S.loc[mic_id, e.id] = stoic
    return S
