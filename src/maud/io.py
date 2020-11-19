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

from typing import Dict, List

import toml

from maud import validation
from maud.data_model import (
    Compartment,
    Drain,
    Enzyme,
    Experiment,
    ExperimentSet,
    KineticModel,
    MaudInput,
    Measurement,
    Metabolite,
    MetaboliteInCompartment,
    Modifier,
    Prior,
    PriorSet,
    Reaction,
    StanCodeSet,
)
from maud.utils import codify


MECHANISM_TO_PARAM_IDS = {
    "uniuni": ["Kcat1", "Kcat2", "Ka"],
    "ordered_unibi": ["Kcat1", "Kcat2", "Ka", "Kp", "Kq", "Kia"],
    "ordered_bibi": ["Kcat1", "Kcat2", "Ka", "Kb", "Kp", "Kq", "Kib", "Kiq"],
    "ping_pong": ["Kcat1", "Kcat2", "Ka", "Kb", "Kp", "Kq", "Kia", "Kib", "Kiq"],
    "ordered_terbi": [
        "Kcat1",
        "Kcat2",
        "Ka",
        "Kb",
        "Kc",
        "Kq",
        "Kia",
        "Kib",
        "Kic",
        "Kiq",
    ],
}


def load_kinetic_model_from_toml(
    parsed_toml: dict, model_id: str = "mi"
) -> KineticModel:
    """Turn the output of toml.load into a KineticModel object.

    :param parsed_toml: Result of running toml.load on a suitable toml file
    :param model_id: String identifying the model.

    """
    compartments = {
        c["id"]: Compartment(id=c["id"], name=c["name"], volume=c["volume"])
        for c in parsed_toml["compartments"]
    }
    metabolites = {
        m["id"]: Metabolite(id=m["id"], name=m["name"])
        for m in parsed_toml["metabolites"]
    }
    mics = {
        f"{m['id']}_{m['compartment']}": MetaboliteInCompartment(
            id=f"{m['id']}_{m['compartment']}",
            metabolite_id=m["id"],
            compartment_id=m["compartment"],
            balanced=m["balanced"],
        )
        for m in parsed_toml["metabolites"]
    }
    reactions = {r["id"]: load_reaction_from_toml(r) for r in parsed_toml["reactions"]}
    if "drains" in parsed_toml.keys():
        drains = {d["id"]: load_drain_from_toml(d) for d in parsed_toml["drains"]}
    else:
        drains = None
    return KineticModel(
        model_id=model_id,
        metabolites=metabolites,
        compartments=compartments,
        mics=mics,
        reactions=reactions,
        drains=drains,
    )


def get_stan_codes(km: KineticModel, experiments: ExperimentSet) -> StanCodeSet:
    """Get the stan codes for a Maud input.

    :param km: KineticModel object
    :param experiments: dictionary mapping experiment ids to Experiment objects
    """
    rxns = km.reactions.values()
    enzyme_ids = [eid for rxn in rxns for eid in rxn.enzymes.keys()]
    mic_codes = codify(km.mics.keys())
    balanced_mic_codes = {
        mic_id: code for mic_id, code in mic_codes.items() if km.mics[mic_id].balanced
    }
    unbalanced_mic_codes = {
        mic_id: code
        for mic_id, code in mic_codes.items()
        if not km.mics[mic_id].balanced
    }
    if km.drains is not None:
        drain_codes = codify(km.drains.keys())
    else:
        drain_codes = {}
    return StanCodeSet(
        metabolite_codes=codify(km.metabolites.keys()),
        mic_codes=mic_codes,
        balanced_mic_codes=balanced_mic_codes,
        unbalanced_mic_codes=unbalanced_mic_codes,
        reaction_codes=codify(km.reactions.keys()),
        experiment_codes=codify([e.id for e in experiments.experiments]),
        enzyme_codes=codify(enzyme_ids),
        drain_codes=drain_codes,
    )


def load_drain_from_toml(toml_drain: dict) -> Drain:
    """Turn a dictionary representing a drain into a Drain object.

    :param toml_drain: Dictionary representing a drain, typically one of
    the values of the 'drains' field in the output of toml.load.

    """

    return Drain(
        id=toml_drain["id"],
        name=toml_drain["name"],
        stoichiometry=toml_drain["stoichiometry"],
    )


def load_reaction_from_toml(toml_reaction: dict) -> Reaction:
    """Turn a dictionary representing a reaction into a Reaction object.

    :param toml_reaction: Dictionary representing a reaction, typically one of
    the values of the 'reactions' field in the output of toml.load.

    """
    enzymes = {}
    subunits = 1
    reversible = (
        toml_reaction["reversible"] if "reversible" in toml_reaction.keys() else None
    )
    is_exchange = (
        toml_reaction["is_exchange"] if "is_exchange" in toml_reaction.keys() else None
    )
    water_stoichiometry = (
        toml_reaction["water_stoichiometry"]
        if "water_stoichiometry" in toml_reaction.keys()
        else 0
    )
    for e in toml_reaction["enzymes"]:
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

        enzymes[e["id"]] = Enzyme(
            id=e["id"],
            name=e["name"],
            reaction_id=toml_reaction["id"],
            modifiers=modifiers,
            subunits=subunits,
        )
    return Reaction(
        id=toml_reaction["id"],
        name=toml_reaction["name"],
        reversible=reversible,
        is_exchange=is_exchange,
        stoichiometry=toml_reaction["stoichiometry"],
        enzymes=enzymes,
        water_stoichiometry=water_stoichiometry,
    )


def get_experiment(raw: Dict) -> Experiment:
    """Extract an Experiment object from a dictionary."""
    out = Experiment(id=raw["id"])
    for target_type in ["metabolite", "reaction", "enzyme"]:
        type_measurements = {}
        for m in raw[target_type + "_measurements"]:
            measurement = Measurement(
                target_id=m["target_id"],
                value=m["value"],
                uncertainty=m["uncertainty"],
                scale="ln",
                target_type=target_type,
            )
            type_measurements.update({m["target_id"]: measurement})
        out.measurements.update({target_type: type_measurements})
    if "knockouts" in raw.keys():
        out.knockouts = raw["knockouts"]
    return out


def extract_priors(list_of_prior_dicts: List[Dict], id_func):
    """Get a list of Prior objects from a list of dictionaries."""
    return [Prior(id_func(p), **p) for p in list_of_prior_dicts]


def load_maud_input_from_toml(filepath: str, id: str = "mi") -> MaudInput:
    """
    Load an MaudInput object from a suitable toml file.

    :param filepath: path to a toml file
    :param id: id for the output object

    """
    parsed_toml = toml.load(filepath)
    kinetic_model = load_kinetic_model_from_toml(parsed_toml, id)
    experiments = ExperimentSet([get_experiment(e) for e in parsed_toml["experiments"]])
    prior_dict = parsed_toml["priors"]
    for k in [
        "inhibition_constants", "tense_dissociation_constants",
        "relaxed_dissociation_constants", "transfer_constants"
    ]:
        if k not in prior_dict.keys():
            prior_dict[k] = {}
    priors = PriorSet(
        km_priors=extract_priors(
            prior_dict["kms"], lambda p: f"km_{p['enzyme_id']}_{p['mic_id']}"
        ),
        kcat_priors=extract_priors(
            prior_dict["kcats"], lambda p: f"kcat_{p['enzyme_id']}"
        ),
        formation_energy_priors=extract_priors(
            prior_dict["formation_energies"],
            lambda p: f"formation_energy_{p['metabolite_id']}",
        ),
        inhibition_constant_priors=extract_priors(
            prior_dict["inhibition_constants"],
            lambda p: f"ki_{p['enzyme_id']}_{p['mic_id']}",
        ),
        relaxed_dissociation_constant_priors=extract_priors(
            prior_dict["relaxed_dissociation_constants"],
            lambda p: f"diss_r_{p['enzyme_id']}_{p['mic_id']}",
        ),
        tense_dissociation_constant_priors=extract_priors(
            prior_dict["tense_dissociation_constants"],
            lambda p: f"diss_r_{p['enzyme_id']}_{p['mic_id']}",
        ),
        transfer_constant_priors=extract_priors(
            prior_dict["transfer_constants"],
            lambda p: f"transfer_constant_{p['enzyme_id']}",
        ),
        unbalanced_metabolite_priors=[
            Prior(
                id=f"{e['id']}_{m['target_id']}",
                location=m["value"],
                scale=m["uncertainty"],
                experiment_id=e["id"],
                mic_id=m["target_id"],
            )
            for e in parsed_toml["experiments"]
            for m in e["metabolite_measurements"]
            if not kinetic_model.mics[m["target_id"]].balanced
        ],
        drain_priors=[
            Prior(
                id=f"{dd['id']}_{e['id']}",
                location=edd["location"],
                scale=edd["scale"],
                drain_id=dd["id"],
                experiment_id=e["id"],
            )
            for dd in parsed_toml["drains"]
            for e in parsed_toml["experiments"]
            for edd in e["drains"]
            if edd["id"] == dd["id"]
        ]
        if "drains" in parsed_toml.keys()
        else [],
    )
    stan_codes = get_stan_codes(kinetic_model, experiments)
    mi = MaudInput(
        experiments=experiments,
        kinetic_model=kinetic_model,
        priors=priors,
        stan_codes=stan_codes,
    )
    validation.validate_maud_input(mi)
    return mi
