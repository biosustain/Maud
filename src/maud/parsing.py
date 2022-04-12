"""Functions for turning raw data into Maud objects."""

from datetime import datetime
from typing import Callable, List, Union

import pandas as pd

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
from maud.data_model.maud_config import MaudConfig, ODEConfig
from maud.data_model.measurement_set import (
    EnzymeKnockout,
    Experiment,
    MeasurementSet,
    MeasurementType,
    PhosphorylationKnockout,
)
from maud.utils import read_with_fallback


def parse_measurement_set(
    raw_measurement_table: pd.DataFrame, raw_biological_config: dict
) -> MeasurementSet:
    """Parse a measurements dataframe.

    :param measurement_table: result of running pd.read_csv on suitable file
    :param raw_biological_config: result of running toml.load on suitable file
    """
    experiments = [
        Experiment(id=e["id"], is_train=e["is_train"], is_test=e["is_test"])
        for e in raw_biological_config["experiment"]
    ]
    y = {
        mt: raw_measurement_table.loc[
            lambda df: df["measurement_type"] == mt.value
        ]
        for mt in MeasurementType
    }
    enz_knockouts = (
        [
            EnzymeKnockout(
                experiment_id=eko["experiment_id"], enzyme_id=eko["enzyme_id"]
            )
            for eko in raw_biological_config["enzyme_knockout"]
        ]
        if "enzyme_knockout" in raw_biological_config.keys()
        else None
    )
    phos_knockouts = (
        [
            PhosphorylationKnockout(
                experiment_id=pko["experiment_id"], enzyme_id=pko["enzyme_id"]
            )
            for pko in raw_biological_config["phos_knockout"]
        ]
        if "phos_knockout" in raw_biological_config.keys()
        else None
    )
    return MeasurementSet(
        yconc=y[MeasurementType.MIC],
        yflux=y[MeasurementType.FLUX],
        yenz=y[MeasurementType.ENZYME],
        enz_knockouts=enz_knockouts,
        phos_knockouts=phos_knockouts,
        experiments=experiments,
    )


def parse_config(raw: dict) -> MaudConfig:
    """Get a MaudConfig object from the result of toml.load.

    :param raw: result of running toml.load on a suitable file
    """
    init_kwargs = raw.copy()
    if "ode_config" in init_kwargs.keys():
        init_kwargs["ode_config"] = ODEConfig(**init_kwargs["ode_config"])
    return MaudConfig(**init_kwargs)


def parse_kinetic_model(raw: dict) -> KineticModel:
    """Turn the output of toml.load into a KineticModel object.

    :param raw: Result of running toml.load on a suitable toml file

    """
    now = datetime.now().strftime("%Y%m%d%H%M%S")
    name = read_with_fallback("name", raw, now)
    compartments = [
        Compartment(c["id"], c["name"], c["volume"])
        for c in raw["compartment"]
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
            mechanism=ReactionMechanism[r["mechanism"].upper()],
            stoichiometry=r["stoichiometry"],
            water_stoichiometry=read_with_fallback(
                "water_stoichiometry", r, 0),
        )
        for r in raw["reaction"]
    ]
    metabolites = parse_metabolites(raw)
    enzymes = parse_enzymes(raw)
    phosphorylations = parse_field_if_available(
        raw,
        "phosphorylation",
        lambda p: Phosphorylation(
            enzyme_id=p["enzyme_id"],
            modification_type=ModificationType[p["modification_type"].upper()],
        ),
    )
    allosteries = parse_field_if_available(
        raw,
        "allostery",
        lambda a: Allostery(
            enzyme_id=a["enzyme_id"],
            metabolite_id=a["metabolite_id"],
            compartment_id=a["compartment_id"],
            modification_type=ModificationType[a["modification_type"].upper()],
        ),
    )
    competitive_inhibitions = parse_field_if_available(
        raw,
        "competitive_inhibition",
        lambda ci: CompetitiveInhibition(
            enzyme_id=ci["enzyme_id"],
            reaction_id=ci["reaction_id"],
            metabolite_id=ci["metabolite_id"],
            compartment_id=ci["compartment_id"],
        ),
    )
    return KineticModel(
        name=name,
        metabolites=metabolites,
        enzymes=enzymes,
        compartments=compartments,
        reactions=reactions,
        mics=mics,
        ers=ers,
        allosteries=allosteries,
        competitive_inhibitions=competitive_inhibitions,
        phosphorylations=phosphorylations,
    )


def parse_metabolites(raw: dict) -> List[Metabolite]:
    """parse explicitly provided metabolites (if any), then infer extra ones.

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
    """parse explicitly provided enzymes (if any), then infer extra ones.

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


def parse_field_if_available(
    raw: dict, field_name: str, parsing_function: Callable
) -> Union[List, None]:
    """Apply a parsing function or return None."""
    if field_name not in raw.keys():
        return None
    else:
        return [parsing_function(x) for x in raw[field_name]]
