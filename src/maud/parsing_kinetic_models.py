"""Functions for parsing kinetic models from raw Maud input data."""

from datetime import datetime
from typing import List

from maud.data_model.kinetic_model import (
    Allostery,
    Compartment,
    CompetitiveInhibition,
    Enzyme,
    EnzymeReaction,
    KineticModel,
    Metabolite,
    MetaboliteInCompartment,
    ModificationType,
    Phosphorylation,
    Reaction,
    ReactionMechanism,
)
from maud.utils import read_with_fallback


def parse_kinetic_model(raw: dict) -> KineticModel:
    """Turn the output of toml.load into a KineticModel object.

    :param raw: Result of running toml.load on a suitable toml file

    """
    now = datetime.now().strftime("%Y%m%d%H%M%S")
    name = read_with_fallback("name", raw, now)
    compartments = [
        Compartment(c["id"], c["name"], c["volume"]) for c in raw["compartment"]
    ]
    mics = [
        MetaboliteInCompartment(
            mic["metabolite_id"], mic["compartment_id"], mic["balanced"]
        )
        for mic in raw["metabolite_in_compartment"]
    ]
    ers = [
        EnzymeReaction(er["enzyme_id"], er["reaction_id"])
        for er in raw["enzyme_reaction"]
    ]
    reactions = [
        Reaction(
            id=r["id"],
            name=read_with_fallback("name", r, None),
            mechanism=ReactionMechanism[r["mechanism"]],
            stoichiometry=r["stoichiometry"],
            water_stoichiometry=read_with_fallback("water_stoichiometry", r, 0),
            transported_charge=read_with_fallback("transported_charge", r, 0),
        )
        for r in raw["reaction"]
    ]
    metabolites = parse_metabolites(raw)
    enzymes = parse_enzymes(raw)
    if "phosphorylation" in raw.keys():
        phosphorylations = [
            Phosphorylation(
                modifying_enzyme_id=p["modifying_enzyme_id"],
                modified_enzyme_id=p["modified_enzyme_id"],
                modification_type=ModificationType[p["modification_type"]],
                name=p["name"] if "name" in p.keys() else None,
            )
            for p in raw["phosphorylation"]
        ]
    else:
        phosphorylations = None
    if "allostery" in raw.keys():
        allosteries = [
            Allostery(
                enzyme_id=a["enzyme_id"],
                metabolite_id=a["metabolite_id"],
                compartment_id=a["compartment_id"],
                modification_type=ModificationType[a["modification_type"]],
            )
            for a in raw["allostery"]
        ]
        allosteric_enzymes = [
            e for e in enzymes if any(a.enzyme_id == e.id for a in allosteries)
        ]
    else:
        allosteries = None
        allosteric_enzymes = None
    if "competitive_inhibition" in raw.keys():
        competitive_inhibitions = [
            CompetitiveInhibition(
                enzyme_id=ci["enzyme_id"],
                reaction_id=ci["reaction_id"],
                metabolite_id=ci["metabolite_id"],
                compartment_id=ci["compartment_id"],
            )
            for ci in raw["competitive_inhibition"]
        ]
    else:
        competitive_inhibitions = None
    return KineticModel(
        name=name,
        metabolites=metabolites,
        enzymes=enzymes,
        compartments=compartments,
        reactions=reactions,
        mics=mics,
        ers=ers,
        allosteries=allosteries,
        allosteric_enzymes=allosteric_enzymes,
        competitive_inhibitions=competitive_inhibitions,
        phosphorylations=phosphorylations,
    )


def parse_metabolites(raw: dict) -> List[Metabolite]:
    """Parse explicitly provided metabolites (if any), then infer extra ones.

    :param raw: output of running toml.load on a suitable file

    """
    metabolites = []
    explicit_met_ids = []
    if "metabolite" in raw.keys():
        metabolites = [
            Metabolite(
                id=m["id"],
                name=read_with_fallback("name", m, None),
                inchi_key=read_with_fallback("inchi_key", m, None),
            )
            for m in raw["metabolite"]
        ]
        explicit_met_ids = [m.id for m in metabolites]
    for mic_met_id in set(
        mic["metabolite_id"] for mic in raw["metabolite_in_compartment"]
    ):
        if mic_met_id not in explicit_met_ids:
            metabolites.append(
                Metabolite(id=mic_met_id, name=None, inchi_key=None)
            )
    return metabolites


def parse_enzymes(raw: dict) -> List[Enzyme]:
    """Parse explicitly provided enzymes (if any), then infer extra ones.

    :param raw: output of running toml.load on a suitable file

    """
    enzymes = []
    explicit_enz_ids = []
    if "enzyme" in raw.keys():
        enzymes = [
            Enzyme(
                id=e["id"],
                name=read_with_fallback("name", e, None),
                subunits=read_with_fallback("subunits", e, 1),
            )
            for e in raw["enzyme"]
        ]
        explicit_enz_ids = [e.id for e in enzymes]
    for er_enz_id in set(er["enzyme_id"] for er in raw["enzyme_reaction"]):
        if er_enz_id not in explicit_enz_ids:
            enzymes.append(Enzyme(id=er_enz_id, name=None, subunits=1))
    return enzymes
