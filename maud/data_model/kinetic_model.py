"""Maud's definition of a kinetic model."""

from enum import Enum
from typing import Dict, List, Optional, Union

import pandas as pd
from pydantic import (
    BaseModel,
    ConfigDict,
    computed_field,
    field_validator,
    model_validator,
)

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


class Metabolite(BaseModel):
    """Maud representation of a metabolite."""

    id: str
    name: Optional[str]
    inchi_key: Optional[str]

    @field_validator("id")
    def id_must_not_contain_seps(cls, v):
        """Check that the id doesn't contain ID_SEPARATOR."""
        assert ID_SEPARATOR not in v
        return v


class Enzyme(BaseModel):
    """Maud representation of an enzyme."""

    id: str
    name: Optional[str]
    subunits: int

    @field_validator("id")
    def id_must_not_contain_seps(cls, v):
        """Check that the id doesn't contain ID_SEPARATOR."""
        assert ID_SEPARATOR not in v
        return v

    @field_validator("subunits")
    def subunits_must_be_positive(cls, v):
        """Check that the subunits attribute is biologically possible."""
        assert v > 0
        return v


class PhosphorylationModifyingEnzyme(BaseModel):
    """Maud representation of a phosphorylation modifying enzyme.

    For example, a phosphatase?

    """

    id: str

    @field_validator("id")
    def id_must_not_contain_seps(cls, v):
        """Check that the id doesn't contain ID_SEPARATOR."""
        assert ID_SEPARATOR not in v
        return v


class Compartment(BaseModel):
    """Maud representation of an intra-cellular compartment.

    For example, cytosol or mitochondria.

    """

    id: str
    name: Optional[str]
    volume: float

    @field_validator("id")
    def id_must_not_contain_seps(cls, v):
        """Check that the id doesn't contain ID_SEPARATOR."""
        assert ID_SEPARATOR not in v
        return v


class Reaction(BaseModel):
    """Maud representation of a chemical reaction."""

    id: str
    name: Optional[str]
    mechanism: ReactionMechanism
    stoichiometry: Dict[str, float]
    water_stoichiometry: float
    transported_charge: float

    @field_validator("id")
    def id_must_not_contain_seps(cls, v):
        """Check that the id doesn't contain ID_SEPARATOR."""
        assert ID_SEPARATOR not in v, "ID must not contain separator"
        return v

    @field_validator("stoichiometry")
    def stoichiometry_must_be_non_zero(cls, v):
        """Check that the stoichiometry is not zero."""
        assert v != 0, "stoichiometry must be non-zero"
        return v


class MetaboliteInCompartment(BaseModel):
    """Maud representation of a metabolite/compartment pair.

    This is needed because metabolites often exist in multiple compartments, and
    the concentration in each one is important.

    A metabolite may also be "balanced" (i.e. where the inflows and outflows are
    accounted for and equal) in one compartment but not in another.

    """

    metabolite_id: str
    compartment_id: str
    balanced: bool

    @computed_field
    def id(self) -> str:
        """Add the id field."""
        return ID_SEPARATOR.join([self.metabolite_id, self.compartment_id])


class EnzymeReaction(BaseModel):
    """Maud representation of an enzyme/reaction pair.

    This is needed because some enzymes catalyse multiple reactions.

    """

    enzyme_id: str
    reaction_id: str

    @computed_field
    def id(self) -> str:
        """Add the id field."""
        return self.enzyme_id + ID_SEPARATOR + self.reaction_id


class Allostery(BaseModel):
    """Maud representation of an allosteric modification."""

    enzyme_id: str
    metabolite_id: str
    compartment_id: str
    modification_type: ModificationType

    @computed_field
    def id(self) -> str:
        """Add the id field."""
        return ID_SEPARATOR.join(
            [
                self.enzyme_id,
                self.metabolite_id,
                self.compartment_id,
                self.modification_type.name,
            ]
        )

    @computed_field
    def mic_id(self) -> str:
        """Add the mic_id field."""
        return ID_SEPARATOR.join([self.metabolite_id, self.compartment_id])


class CompetitiveInhibition(BaseModel):
    """Maud representation of a competitive inhibition."""

    enzyme_id: str
    reaction_id: str
    metabolite_id: str
    compartment_id: str

    @computed_field
    def id(self) -> str:
        """Add the id field."""
        return ID_SEPARATOR.join(
            [
                self.enzyme_id,
                self.reaction_id,
                self.metabolite_id,
                self.compartment_id,
            ]
        )

    @computed_field
    def er_id(self) -> str:
        """Add the er_id field."""
        return ID_SEPARATOR.join([self.enzyme_id, self.reaction_id])

    @computed_field
    def mic_id(self) -> str:
        """Add the mic_id field."""
        return ID_SEPARATOR.join([self.metabolite_id, self.compartment_id])


class Phosphorylation(BaseModel):
    """Maud representation of a phosphorylation modification."""

    name: Optional[str]
    modifying_enzyme_id: str
    modified_enzyme_id: str
    modification_type: ModificationType

    @computed_field
    def id(self) -> str:
        """Add the id field."""
        return ID_SEPARATOR.join(
            [
                self.modifying_enzyme_id,
                self.modified_enzyme_id,
                self.modification_type.name,
            ]
        )


class KineticModel(BaseModel):
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
    model_config: ConfigDict = {"arbitrary_types_allowed": True}

    @computed_field
    def drains(self) -> List[Reaction]:
        """Add the drains field."""
        return [
            r for r in self.reactions if r.mechanism == ReactionMechanism.drain
        ]

    @computed_field
    def edges(self) -> List[Union[Reaction, EnzymeReaction]]:
        """Add the edges field."""
        return self.drains + self.ers

    @computed_field
    def stoichiometric_matrix(self) -> pd.DataFrame:
        """Add the stoichiometric_matrix field."""
        return get_stoichiometric_matrix(self.edges, self.mics, self.reactions)

    @computed_field
    def phosphorylation_modifying_enzymes(
        self,
    ) -> Optional[List[PhosphorylationModifyingEnzyme]]:
        """Add the phosphorylation_modifying_enzymes field."""
        return (
            [
                PhosphorylationModifyingEnzyme(pme_id)
                for pme_id in list(
                    set(p.modifying_enzyme_id for p in self.phosphorylations)
                )
            ]
            if self.phosphorylations is not None
            else None
        )

    @field_validator("metabolites")
    def metabolite_ids_must_be_unique(cls, v):
        """Make sure there aren't any duplicated metabolite ids."""
        met_ids = [m.id for m in v]
        assert len(met_ids) == len(set(met_ids))
        return v

    @field_validator("enzymes")
    def enzyme_ids_must_be_unique(cls, v):
        """Make sure there aren't any duplicated enzyme ids."""
        met_ids = [m.id for m in v]
        assert len(met_ids) == len(set(met_ids))
        return v

    @field_validator("compartments")
    def compartment_ids_must_be_unique(cls, v):
        """Make sure there aren't any duplicated compartment ids."""
        met_ids = [m.id for m in v]
        assert len(met_ids) == len(set(met_ids))
        return v

    @field_validator("reactions")
    def reaction_ids_must_be_unique(cls, v):
        """Make sure there aren't any duplicated reaction ids."""
        rxn_ids = [r.id for r in v]
        assert len(rxn_ids) == len(set(rxn_ids))
        return v

    @model_validator(mode="after")
    def stoic_keys_must_be_mic_ids(self) -> "KineticModel":
        """Make sure reaction stoichiometries have existent mic ids."""
        mic_ids = [mic.id for mic in self.mics]
        for r in self.reactions:
            for stoich_mic_id in r.stoichiometry.keys():
                assert (
                    stoich_mic_id in mic_ids
                ), f"{r.id} has stoichiometry for bad mic_id {stoich_mic_id}"
        return self

    @model_validator(mode="after")
    def mic_references_must_exist(self) -> "KineticModel":
        """Make sure mics have existent metabolite and compartment ids."""
        metabolite_ids = [m.id for m in self.metabolites]
        compartment_ids = [c.id for c in self.compartments]
        for mic in self.mics:
            assert (
                mic.metabolite_id in metabolite_ids
            ), f"{mic.id} has bad metabolite_id."
            assert (
                mic.compartment_id in compartment_ids
            ), f"{mic.id} has bad compartment_id."
        return self

    @model_validator(mode="after")
    def er_references_must_exist(self) -> "KineticModel":
        """Make sure ers have existent enzyme and reaction ids."""
        enzyme_ids = [e.id for e in self.enzymes]
        reaction_ids = [r.id for r in self.reactions]
        for er in self.ers:
            assert er.enzyme_id in enzyme_ids, f"{er.id} has bad enzyme_id"
            assert (
                er.reaction_id in reaction_ids
            ), f"{er.id} has bad reaction_id"
        return self

    @model_validator(mode="after")
    def allostery_references_must_exist(self) -> "KineticModel":
        """Make sure allosteries' external ids exist."""
        if self.allosteries is None:
            return self
        enzyme_ids = [e.id for e in self.enzymes]
        mic_ids = [mic.id for mic in self.mics]
        for allostery in self.allosteries:
            assert (
                allostery.enzyme_id in enzyme_ids
            ), f"{allostery.id} has bad enzyme_id"
            assert allostery.mic_id in mic_ids, f"{allostery.id} has bad mic_id"
        return self

    @model_validator(mode="after")
    def ci_references_must_exist(self) -> "KineticModel":
        """Make sure competitive inhibitions' external ids exist."""
        if self.competitive_inhibitions is None:
            return self
        er_ids = [er.id for er in self.ers]
        mic_ids = [mic.id for mic in self.mics]
        for ci in self.competitive_inhibitions:
            assert ci.er_id in er_ids, f"{ci.id} has bad er_id"
            assert ci.mic_id in mic_ids, f"{ci.mic_id} has bad mic_id"
        return self

    @model_validator(mode="after")
    def phosphorylation_references_must_exist(self):
        """Make sure phosphorylations' external ids exist."""
        if self.phosphorylations is None:
            return self
        enzyme_ids = [e.id for e in self.enzymes]
        for p in self.phosphorylations:
            assert (
                p.modified_enzyme_id in enzyme_ids
            ), f"{p.id} has bad enzyme_id"
        return self


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
