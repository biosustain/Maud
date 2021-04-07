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

"""Functions for loading MaudInput objects.

(and at some point in the future also saving MaudOutput objects)

"""

import os
from typing import List, Tuple

import pandas as pd
import toml

from maud import validation
from maud.data_model import (
    Compartment,
    Drain,
    Enzyme,
    KineticModel,
    Knockout,
    MaudConfig,
    MaudInput,
    Measurement,
    Metabolite,
    MetaboliteInCompartment,
    Modifier,
    Phosphorylation,
    Prior,
    PriorSet,
    Reaction,
    StanCoordSet,
)
from maud.utils import codify


def load_maud_input_from_toml(data_path: str) -> MaudInput:
    """
    Load an MaudInput object from a data path.

    :param filepath: path to directory containing input toml file

    """
    config = parse_config(toml.load(os.path.join(data_path, "config.toml")))
    kinetic_model_path = os.path.join(data_path, config.kinetic_model_file)
    experiments_path = os.path.join(data_path, config.experiments_file)
    priors_path = os.path.join(data_path, config.priors_file)
    kinetic_model = parse_toml_kinetic_model(toml.load(kinetic_model_path))
    measurements, knockouts = parse_measurements(pd.read_csv(experiments_path))
    stan_coords = get_stan_coords(kinetic_model, measurements)
    prior_set = parse_prior_set_df(pd.read_csv(priors_path), stan_coords)
    mi = MaudInput(
        config=config,
        kinetic_model=kinetic_model,
        priors=prior_set,
        stan_coords=stan_coords,
        measurements=measurements,
        knockouts=knockouts,
    )
    validation.validate_maud_input(mi)
    return mi


def parse_toml_kinetic_model(raw: dict) -> KineticModel:
    """Turn the output of toml.load into a KineticModel object.

    :param raw: Result of running toml.load on a suitable toml file

    """
    model_id = raw["model_id"] if "model_id" in raw.keys() else "mi"
    compartments = [
        Compartment(id=c["id"], name=c["name"], volume=c["volume"])
        for c in raw["compartments"]
    ]
    raw_mets = {m["id"]: m for m in raw["metabolites"]}
    metabolites = [Metabolite(id=m["id"], name=m["name"]) for m in raw_mets.values()]
    mics = [
        MetaboliteInCompartment(
            id=f"{m['id']}_{m['compartment']}",
            metabolite_id=m["id"],
            compartment_id=m["compartment"],
            balanced=m["balanced"],
        )
        for m in raw["metabolites"]
    ]
    reactions = [parse_toml_reaction(r) for r in raw["reactions"]]
    drains = (
        [parse_toml_drain(d) for d in raw["drains"]] if "drains" in raw.keys() else []
    )
    phosphorylation = (
        [parse_toml_phosphorylation(d) for d in raw["phosphorylation"]]
        if "phosphorylation" in raw.keys()
        else []
    )
    return KineticModel(
        model_id=model_id,
        metabolites=metabolites,
        compartments=compartments,
        mics=mics,
        reactions=reactions,
        drains=drains,
        phosphorylation=phosphorylation,
    )


def get_stan_coords(km: KineticModel, ms: List[Measurement]) -> StanCoordSet:
    """Get the stan codes for a Maud input.

    :param km: KineticModel object
    :param ms: MeasurementSet object
    """
    kms = [
        f"{e.id}_{m}" for r in km.reactions for e in r.enzymes for m in r.stoichiometry
    ]
    experiments = sorted(
        list(set(m.experiment_id for m in ms))
    )  # sorted because set has non-deterministic order
    mics = [m.id for m in km.mics]
    reactions = [r.id for r in km.reactions]
    enzymes = [e.id for r in km.reactions for e in r.enzymes]
    yconc_exps, yflux_exps, yenz_exps = (
        [m.experiment_id for m in ms if m.target_type == t]
        for t in ["mic", "flux", "enzyme"]
    )
    yconc_mics, yflux_rxns, yenz_enzs = (
        [m.target_id for m in ms if m.target_type == t]
        for t in ["mic", "flux", "enzyme"]
    )
    ci_mics, ai_mics, aa_mics = (
        [
            m.mic_id
            for r in km.reactions
            for e in r.enzymes
            for m in e.modifiers[modifier_type]
        ]
        for modifier_type in [
            "competitive_inhibitor",
            "allosteric_inhibitor",
            "allosteric_activator",
        ]
    )
    return StanCoordSet(
        metabolites=[m.id for m in km.metabolites],
        mics=mics,
        kms=kms,
        balanced_mics=[m.id for m in km.mics if m.balanced],
        unbalanced_mics=[m.id for m in km.mics if not m.balanced],
        reactions=reactions,
        experiments=experiments,
        enzymes=enzymes,
        phos_enzs=[e.id for e in km.phosphorylation],
        drains=[d.id for d in km.drains],
        yconc_exps=yconc_exps,
        yconc_mics=yconc_mics,
        yflux_exps=yflux_exps,
        yflux_rxns=yflux_rxns,
        yenz_exps=yenz_exps,
        yenz_enzs=yenz_enzs,
        ci_mics=ci_mics,
        ai_mics=ai_mics,
        aa_mics=aa_mics,
    )


def parse_toml_drain(raw: dict) -> Drain:
    """Turn a dictionary representing a drain into a Drain object.

    :param raw: Dictionary representing a drain, typically one of
    the values of the 'drains' field in the output of toml.load.

    """

    return Drain(id=raw["id"], name=raw["name"], stoichiometry=raw["stoichiometry"])


def parse_toml_phosphorylation(raw: dict) -> Phosphorylation:
    """Turn a dictionary-format phosphorylation into a Maud object.

    :param toml_phosphorylation: Dictionary representing a phosphorylation
    enzyme, typically one of the values of the 'phosphorylation' field in the
    output of toml.load.

    """

    return Phosphorylation(
        id=raw["id"],
        name=raw["name"],
        activating=raw["activating"] if "activating" in raw.keys() else None,
        inhibiting=raw["inhibiting"] if "inhibiting" in raw.keys() else None,
        enzyme_id=raw["enzyme_id"],
    )


def parse_toml_reaction(raw: dict) -> Reaction:
    """Turn a dictionary representing a reaction into a Reaction object.

    :param raw: Dictionary representing a reaction, typically one of
    the values of the 'reactions' field in the output of toml.load.

    """
    enzymes = []
    subunits = 1
    water_stoichiometry = (
        raw["water_stoichiometry"] if "water_stoichiometry" in raw.keys() else 0
    )
    for e in raw["enzymes"]:
        modifiers = {
            "competitive_inhibitor": [],
            "allosteric_inhibitor": [],
            "allosteric_activator": [],
        }
        if "modifiers" in e.keys():
            for modifier_dict in e["modifiers"]:
                modifier_type = modifier_dict["modifier_type"]
                modifiers[modifier_type].append(
                    Modifier(
                        mic_id=modifier_dict["mic_id"],
                        enzyme_id=e["id"],
                        modifier_type=modifier_type,
                    )
                )

        if "subunits" in e.keys():
            subunits = e["subunits"]

        enzymes.append(
            Enzyme(
                id=e["id"],
                name=e["name"],
                reaction_id=raw["id"],
                modifiers=modifiers,
                subunits=subunits,
                water_stoichiometry=water_stoichiometry,
            )
        )
    return Reaction(
        id=raw["id"],
        name=raw["name"],
        stoichiometry=raw["stoichiometry"],
        enzymes=enzymes,
    )


def parse_measurements(
    raw: pd.DataFrame,
) -> Tuple[List[Measurement], List[Knockout]]:
    """Parse a measurements dataframe.

    :param raw: result of running pd.read_csv on a suitable file
    """
    measurements = []
    knockouts = []
    for _, row in raw.iterrows():
        if row["measurement_type"] in ["knockout_enz", "knockout_phos"]:
            knockouts.append(
                Knockout(
                    experiment_id=row["experiment_id"],
                    target_id=row["target_id"],
                    knockout_type=row["measurement_type"],
                )
            )
        else:
            measurements.append(
                Measurement(
                    target_id=row["target_id"],
                    experiment_id=row["experiment_id"],
                    target_type=row["measurement_type"],
                    value=row["measurement"],
                    error=row["error_scale"],
                )
            )
    return measurements, knockouts


def parse_prior_set_df(raw: pd.DataFrame, cs: StanCoordSet) -> PriorSet:
    """Get a PriorSet object from a dataframe.

    :param raw: result of running pd.read_csv on a suitable file
    """
    priors = {
        "kcat_priors": [],
        "km_priors": [],
        "formation_energy_priors": [],
        "unbalanced_metabolite_priors": [],
        "inhibition_constant_priors": [],
        "tense_dissociation_constant_priors": [],
        "relaxed_dissociation_constant_priors": [],
        "transfer_constant_priors": [],
        "drain_priors": [],
        "enzyme_concentration_priors": [],
        "phos_kcat_priors": [],
        "phos_enz_concentration_priors": [],
    }
    order_keys = {
        # ensure that priors are in the right order wrt the stan codes
        "kcat_priors": lambda p: codify(cs.enzymes)[p.enzyme_id],
        "km_priors": lambda p: codify(cs.kms)[f"{p.enzyme_id}_{p.mic_id}"],
        "formation_energy_priors": lambda p: codify(cs.metabolites)[p.metabolite_id],
        "unbalanced_metabolite_priors": lambda p: (
            codify(cs.experiments)[p.experiment_id],
            codify(cs.mics)[p.mic_id],
        ),
        "inhibition_constant_priors": lambda p: (
            codify(cs.enzymes)[p.enzyme_id],
            codify(cs.ci_mics)[p.mic_id],
        ),
        "tense_dissociation_constant_priors": lambda p: (
            codify(cs.enzymes)[p.enzyme_id],
            codify(cs.ai_mics)[p.mic_id],
        ),
        "relaxed_dissociation_constant_priors": lambda p: (
            codify(cs.enzymes)[p.enzyme_id],
            codify(cs.aa_mics)[p.mic_id],
        ),
        "transfer_constant_priors": lambda p: codify(cs.enzymes)[p.enzyme_id],
        "drain_priors": lambda p: codify(cs.drains)[p.drain_id],
        "enzyme_concentration_priors": lambda p: codify(cs.enzymes)[p.enzyme_id],
        "phos_kcat_priors": lambda p: codify(cs.phos_enzs)[p.phos_enz_id],
        "phos_enz_concentration_priors": lambda p: codify(cs.phos_enzs)[p.phos_enz_id],
    }
    negative_param_types = ["formation_energy", "drain"]
    for _, row in raw.iterrows():
        id = "_".join(row.loc[lambda s: [isinstance(x, str) for x in s]].values)
        parameter_type = row["parameter_type"]
        prior_dict = row.dropna().drop("parameter_type").to_dict()
        is_non_negative = parameter_type not in negative_param_types
        priors[parameter_type + "_priors"].append(
            Prior(id=id, is_non_negative=is_non_negative, **prior_dict)
        )
    for k, v in priors.items():
        priors[k] = sorted(v, key=order_keys[k])
    return PriorSet(**priors)


def parse_config(raw):
    """Get a MaudConfig object from the result of toml.load.

    :param raw: result of running toml.load on a suitable file
    """
    return MaudConfig(
        name=raw["name"],
        kinetic_model_file=raw["kinetic_model"],
        priors_file=raw["priors"],
        experiments_file=raw["experiments"],
        likelihood=raw["likelihood"],
        ode_config=raw["ode_config"],
        cmdstanpy_config=raw["cmdstanpy_config"],
    )
