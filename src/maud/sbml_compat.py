# Copyright (C) 2021 Moritz E. Beber.
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


"""Provide compatibility with SBML."""


from __future__ import annotations

from pathlib import Path
from typing import List, Dict

import libsbml

from .data_model import (
    Compartment,
    Metabolite,
    MetaboliteInCompartment,
    Enzyme,
    Reaction,
    KineticModel,
)


class MaudModelBuilder:
    """Define a builder for the kinetic model aggregate."""

    def __init__(self, **kwargs):
        """"""
        super().__init__(**kwargs)
        self._compartments: List[Compartment] = []
        self._metabolites: List[Metabolite] = []
        self._reactions: List[Reaction] = []
        self._mics: List[MetaboliteInCompartment] = []

    def build_model(self, model_id: str) -> KineticModel:
        """"""
        return KineticModel(
            model_id=model_id,
            metabolites=self._metabolites,
            reactions=self._reactions,
            compartments=self._compartments,
            mics=self._mics,
        )

    def build_compartment(self, id: str, name: str, volume: float = 1.0) -> None:
        """"""
        self._compartments.append(Compartment(id=id, name=name, volume=volume))

    def build_metabolite(
        self, id: str, name: str, compartment_id: str, balanced: bool
    ) -> None:
        """"""
        self._metabolites.append(Metabolite(id=id, name=name))
        self._mics.append(
            MetaboliteInCompartment(
                id=f"{id}_{compartment_id}",
                metabolite_id=id,
                compartment_id=compartment_id,
                balanced=balanced,
            )
        )

    def build_reaction(
        self,
        id: str,
        name: str,
        mechanism: str,
        stoichiometry: Dict[str, float],
        enzymes: List[Enzyme],
    ) -> None:
        """"""
        self._reactions.append(
            Reaction(
                id=id,
                name=name,
                reaction_mechanism=mechanism,
                stoichiometry=stoichiometry,
                enzymes=enzymes,
            )
        )


class SBMLModelParser:
    """Define an SBML model parser which hands off information to a builder instance."""

    def __init__(self, *, document: libsbml.SBMLDocument, **kwargs) -> None:
        """"""
        super().__init__(**kwargs)
        self._doc = document
        self._builder = None

    @classmethod
    def from_file(cls, path: Path) -> SBMLModelParser:
        """"""
        return cls(document=libsbml.readSBMLFromFile(str(path)))

    def parse(self, builder: MaudModelBuilder) -> KineticModel:
        """"""
        self._builder = builder
        model: libsbml.Model = self._doc.getModel()
        self._parse_compartments(model)
        self._parse_metabolites(model)
        self._parse_reactions(model)
        return self._builder.build_model(model.getIdAttribute())

    def _parse_compartments(self, model: libsbml.Model) -> None:
        """"""
        for compartment in model.getListOfCompartments():
            self._builder.build_compartment(
                compartment.getIdAttribute(),
                compartment.getName(),
                compartment.getVolume(),
            )

    def _parse_metabolites(self, model: libsbml.Model) -> None:
        """"""
        for species in model.getListOfSpecies():
            self._builder.build_metabolite(
                species.getIdAttribute(),
                species.getName(),
                species.getCompartment(),
                not species.getBoundaryCondition(),
            )

    def _parse_reactions(self, model: libsbml.Model) -> None:
        """"""
        for reaction in model.getListOfReactions():
            stoichiometry = {}
            for reactant in reaction.getListOfReactants():
                stoichiometry[reactant.getSpecies()] = -float(
                    reactant.getStoichiometry()
                )
            for reactant in reaction.getListOfProducts():
                stoichiometry[reactant.getSpecies()] = float(
                    reactant.getStoichiometry()
                )
            # FIXME: The reaction mechanism could be parsed from an SBO term, e.g.,
            #  https://www.ebi.ac.uk/sbo/main/SBO:0000239 (might have to be added to the
            #  SBO) or a custom annotation. Please ask Matthias König or Andreas Dräger
            #  for advice. Or post on the SBML mailing list.
            enzymes = []
            # FIXME: I must admit, I'm not sure what the best SBML type is to hold the
            #  Maud-specific enzyme information. Maybe you can get some inspiration from
            #  looking at a few kinetic models defined in SBML or the tellurium format
            #  https://tellurium.readthedocs.io/.
            self._builder.build_reaction(
                reaction.getIdAttribute(),
                reaction.getName(),
                "reversible_modular_rate_law",
                stoichiometry,
                enzymes,
            )
