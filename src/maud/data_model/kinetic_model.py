"""Maud's definition of a kinetic model."""

from enum import Enum
from typing import Dict, List, Optional, Union

import pandas as pd
from pydantic import Field, root_validator, validator
from pydantic.dataclasses import dataclass

from maud.data_model.hardcoding import ID_SEPARATOR


class ReactionMechanism(int, Enum):
    """Possible reaction mechanisms."""

    reversible_michaelis_menten = 1
    irreversible_michaelis_menten = 2
    drain = 3


class ModificationType(int, Enum):
    """Possible modification types."""

    activation = 1
    inhibition = 2


class KMConfig:
    """Config allowing the KineticModel class to contain pandas objects."""

    arbitrary_types_allowed = True


@dataclass
class Metabolite:
    """Maud representation of a metabolite."""

    id: str
    name: Optional[str]
    inchi_key: Optional[str]

    @validator("id")
    def id_must_not_contain_seps(cls, v):
        """Check that the id doesn't contain ID_SEPARATOR."""
        assert ID_SEPARATOR not in v
        return v


@dataclass
class Enzyme:
    """Maud representation of an enzyme."""

    id: str
    name: Optional[str]
    subunits: int

    @validator("id")
    def id_must_not_contain_seps(cls, v):
        """Check that the id doesn't contain ID_SEPARATOR."""
        assert ID_SEPARATOR not in v
        return v

    @validator("subunits")
    def subunits_must_be_positive(cls, v):
        """Check that the subunits attribute is biologically possible."""
        assert v > 0
        return v


@dataclass
class PhosphorylationModifyingEnzyme:
    """Maud representation of a phosphorylation modifying enzyme.

    For example, a phosphatase?

    """

    id: str

    @validator("id")
    def id_must_not_contain_seps(cls, v):
        """Check that the id doesn't contain ID_SEPARATOR."""
        assert ID_SEPARATOR not in v
        return v


@dataclass
class Compartment:
    """Maud representation of an intra-cellular compartment.

    For example, cytosol or mitochondria.

    """

    id: str
    name: Optional[str]
    volume: float

    @validator("id")
    def id_must_not_contain_seps(cls, v):
        """Check that the id doesn't contain ID_SEPARATOR."""
        assert ID_SEPARATOR not in v
        return v


@dataclass
class Reaction:
    """Maud representation of a chemical reaction."""

    id: str
    name: Optional[str]
    mechanism: ReactionMechanism
    stoichiometry: Dict[str, float]
    water_stoichiometry: float
    transported_charge: float

    @validator("id")
    def id_must_not_contain_seps(cls, v):
        """Check that the id doesn't contain ID_SEPARATOR."""
        assert ID_SEPARATOR not in v, "ID must not contain separator"
        return v

    @validator("stoichiometry")
    def stoichiometry_must_be_non_zero(cls, v):
        """Check that the stoichiometry is not zero."""
        assert v != 0, "stoichiometry must be non-zero"
        return v


@dataclass
class MetaboliteInCompartment:
    """Maud representation of a metabolite/compartment pair.

    This is needed because metabolites often exist in multiple compartments, and
    the concentration in each one is important.

    A metabolite may also be "balanced" (i.e. where the inflows and outflows are
    accounted for and equal) in one compartment but not in another.

    """

    metabolite_id: str
    compartment_id: str
    balanced: bool
    id: str = Field(init=False, exclude=True)

    def __post_init__(self):
        """Add the id field."""
        self.id = ID_SEPARATOR.join([self.metabolite_id, self.compartment_id])


@dataclass
class EnzymeReaction:
    """Maud representation of an enzyme/reaction pair.

    This is needed because some enzymes catalyse multiple reactions.

    """

    enzyme_id: str
    reaction_id: str
    id: str = Field(init=False, exclude=True)

    def __post_init__(self):
        """Add the id field."""
        self.id = self.enzyme_id + ID_SEPARATOR + self.reaction_id


@dataclass
class Allostery:
    """Maud representation of an allosteric modification."""

    enzyme_id: str
    metabolite_id: str
    compartment_id: str
    modification_type: ModificationType
    id: str = Field(init=False, exclude=True)

    def __post_init__(self):
        """Add the id and mic_id fields."""
        self.id = ID_SEPARATOR.join(
            [
                self.enzyme_id,
                self.metabolite_id,
                self.compartment_id,
                self.modification_type.name,
            ]
        )
        self.mic_id = ID_SEPARATOR.join(
            [self.metabolite_id, self.compartment_id]
        )


@dataclass
class CompetitiveInhibition:
    """Maud representation of a competitive inhibition."""

    enzyme_id: str
    reaction_id: str
    metabolite_id: str
    compartment_id: str
    id: str = Field(init=False, exclude=True)

    def __post_init__(self):
        """Add the id, er_id and mic_id fields."""
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
    """Maud representation of a phosphorylation modification."""

    name: Optional[str]
    modifying_enzyme_id: str
    modified_enzyme_id: str
    modification_type: ModificationType
    id: str = Field(init=False, exclude=True)

    def __post_init__(self):
        """Add the id field."""
        self.id = ID_SEPARATOR.join(
            [
                self.modifying_enzyme_id,
                self.modified_enzyme_id,
                self.modification_type.name,
            ]
        )


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
    allosteric_enzymes: Optional[List[Enzyme]]
    competitive_inhibitions: Optional[List[CompetitiveInhibition]]
    phosphorylations: Optional[List[Phosphorylation]]
    drains: List[Reaction] = Field(init=False, exclude=True)
    edges: List[Union[Reaction, EnzymeReaction]] = Field(
        init=False, exclude=True
    )
    stoichiometric_matrix: pd.DataFrame = Field(init=False, exclude=True)
    phosphorylation_modifying_enzymes: Optional[
        List[PhosphorylationModifyingEnzyme]
    ] = Field(init=False, exclude=True)

    def __post_init__(self):
        """Add drains, edges and stoichiometric matrix."""
        self.drains = [
            r for r in self.reactions if r.mechanism == ReactionMechanism.drain
        ]
        self.edges = self.drains + self.ers
        self.stoichiometric_matrix = get_stoichiometric_matrix(
            self.edges, self.mics, self.reactions
        )

        self.phosphorylation_modifying_enzymes = (
            [
                PhosphorylationModifyingEnzyme(pme_id)
                for pme_id in list(
                    set(p.modifying_enzyme_id for p in self.phosphorylations)
                )
            ]
            if self.phosphorylations is not None
            else None
        )

    @validator("metabolites")
    def metabolite_ids_must_be_unique(cls, v):
        """Make sure there aren't any duplicated metabolite ids."""
        met_ids = [m.id for m in v]
        assert len(met_ids) == len(set(met_ids))
        return v

    @validator("enzymes")
    def enzyme_ids_must_be_unique(cls, v):
        """Make sure there aren't any duplicated enzyme ids."""
        met_ids = [m.id for m in v]
        assert len(met_ids) == len(set(met_ids))
        return v

    @validator("compartments")
    def compartment_ids_must_be_unique(cls, v):
        """Make sure there aren't any duplicated compartment ids."""
        met_ids = [m.id for m in v]
        assert len(met_ids) == len(set(met_ids))
        return v

    @validator("reactions")
    def reaction_ids_must_be_unique(cls, v):
        """Make sure there aren't any duplicated reaction ids."""
        rxn_ids = [r.id for r in v]
        assert len(rxn_ids) == len(set(rxn_ids))
        return v

    @root_validator(pre=False, skip_on_failure=True)
    def stoic_keys_must_be_mic_ids(cls, values):
        """Make sure reaction stoichiometries have existent mic ids."""
        mic_ids = [mic.id for mic in values["mics"]]
        for r in values["reactions"]:
            for stoich_mic_id in r.stoichiometry.keys():
                assert (
                    stoich_mic_id in mic_ids
                ), f"{r.id} has stoichiometry for bad mic_id {stoich_mic_id}"
        return values

    @root_validator(pre=False, skip_on_failure=True)
    def mic_references_must_exist(cls, values):
        """Make sure mics have existent metabolite and compartment ids."""
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

    @root_validator(pre=False, skip_on_failure=True)
    def er_references_must_exist(cls, values):
        """Make sure ers have existent enzyme and reaction ids."""
        enzyme_ids = [e.id for e in values["enzymes"]]
        reaction_ids = [r.id for r in values["reactions"]]
        for er in values["ers"]:
            assert er.enzyme_id in enzyme_ids, f"{er.id} has bad enzyme_id"
            assert (
                er.reaction_id in reaction_ids
            ), f"{er.id} has bad reaction_id"
        return values

    @root_validator(pre=False, skip_on_failure=True)
    def allostery_references_must_exist(cls, values):
        """Make sure allosteries' external ids exist."""
        if values["allosteries"] is None:
            return values
        enzyme_ids = [e.id for e in values["enzymes"]]
        mic_ids = [mic.id for mic in values["mics"]]
        for allostery in values["allosteries"]:
            assert (
                allostery.enzyme_id in enzyme_ids
            ), f"{allostery.id} has bad enzyme_id"
            assert allostery.mic_id in mic_ids, f"{allostery.id} has bad mic_id"
        return values

    @root_validator(pre=False, skip_on_failure=True)
    def ci_references_must_exist(cls, values):
        """Make sure competitive inhibitions' external ids exist."""
        if values["competitive_inhibitions"] is None:
            return values
        er_ids = [er.id for er in values["ers"]]
        mic_ids = [mic.id for mic in values["mics"]]
        for ci in values["competitive_inhibitions"]:
            assert ci.er_id in er_ids, f"{ci.id} has bad er_id"
            assert ci.mic_id in mic_ids, f"{ci.mic_id} has bad mic_id"
        return values

    @root_validator(pre=False, skip_on_failure=True)
    def phosphorylation_references_must_exist(cls, values):
        """Make sure phosphorylations' external ids exist."""
        if values["phosphorylations"] is None:
            return values
        enzyme_ids = [e.id for e in values["enzymes"]]
        for p in values["phosphorylations"]:
            assert (
                p.modified_enzyme_id in enzyme_ids
            ), f"{p.id} has bad enzyme_id"
        return values


def get_stoichiometric_matrix(
    edges: List[Union[Reaction, EnzymeReaction]],
    mics: List[MetaboliteInCompartment],
    rxns: List[Reaction],
) -> pd.DataFrame:
    """Get a stoichiometric matrix from lists of edges, mics and reactions."""
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
