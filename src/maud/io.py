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

from typing import Dict

import toml

from maud import validation
from maud.data_model import (
    Compartment,
    Enzyme,
    Experiment,
    KineticModel,
    MaudInput,
    Measurement,
    Metabolite,
    MetaboliteInCompartment,
    Modifier,
    Prior,
    Reaction,
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
    return KineticModel(
        model_id=model_id,
        metabolites=metabolites,
        compartments=compartments,
        mics=mics,
        reactions=reactions,
    )


def get_stan_codes(km: KineticModel, experiments) -> Dict[str, Dict[str, int]]:
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
    return {
        "metabolite": codify(km.metabolites.keys()),
        "metabolite_in_compartment": mic_codes,
        "balanced_mic": balanced_mic_codes,
        "unbalanced_mic": unbalanced_mic_codes,
        "reaction": codify(km.reactions.keys()),
        "experiment": codify(experiments.keys()),
        "enzyme": codify(enzyme_ids),
    }


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


def load_maud_input_from_toml(filepath: str, id: str = "mi") -> MaudInput:
    """
    Load an MaudInput object from a suitable toml file.

    :param filepath: path to a toml file
    :param id: id for the output object

    """
    parsed_toml = toml.load(filepath)
    kinetic_model = load_kinetic_model_from_toml(parsed_toml, id)
    experiments = {}
    for e in parsed_toml["experiments"]:
        experiment = Experiment(id=e["id"])
        for target_type in ["metabolite", "reaction", "enzyme"]:
            type_measurements = {}
            for m in e[target_type + "_measurements"]:
                measurement = Measurement(
                    target_id=m["target_id"],
                    value=m["value"],
                    uncertainty=m["uncertainty"],
                    scale="ln",
                    target_type=target_type,
                )
                type_measurements.update({m["target_id"]: measurement})
            experiment.measurements.update({target_type: type_measurements})
        if "knockouts" in e.keys():
            experiment.knockouts = e["knockouts"]
        experiments.update({experiment.id: experiment})
    priors = {
        "kcats": [],
        "kms": [],
        "enzyme_concentrations": [],
        "formation_energies": [],
        "unbalanced_metabolites": [],
        "inhibition_constants": [],
        "tense_dissociation_constants": [],
        "relaxed_dissociation_constants": [],
        "transfer_constants": [],
    }
    for prior in parsed_toml["priors"]["formation_energies"]:
        metabolite_id = prior["metabolite_id"]
        priors["formation_energies"].append(
            Prior(
                id=f"formation_energy_{metabolite_id}",
                location=prior["location"],
                scale=prior["scale"],
                metabolite_id=metabolite_id,
            )
        )
    for prior in parsed_toml["priors"]["kms"]:
        enzyme_id = prior["enzyme_id"]
        mic_id = prior["mic_id"]
        priors["kms"].append(
            Prior(
                id=f"km_{enzyme_id}_{mic_id}",
                location=prior["location"],
                scale=prior["scale"],
                mic_id=mic_id,
                enzyme_id=enzyme_id,
            )
        )
    for prior in parsed_toml["priors"]["kcats"]:
        enzyme_id = prior["enzyme_id"]
        priors["kcats"].append(
            Prior(
                id=f"kcat_{enzyme_id}",
                location=prior["location"],
                scale=prior["scale"],
                enzyme_id=enzyme_id,
            )
        )
    for prior_type, prefix in zip(
        [
            "inhibition_constants",
            "relaxed_dissociation_constants",
            "tense_dissociation_constants",
        ],
        ["ki", "diss_r", "diss_t"],
    ):
        if prior_type in parsed_toml["priors"].keys():
            for prior in parsed_toml["priors"][prior_type]:
                enzyme_id = prior["enzyme_id"]
                mic_id = prior["mic_id"]
                priors[prior_type].append(
                    Prior(
                        id=f"{prefix}_{enzyme_id}_{mic_id}",
                        location=prior["location"],
                        scale=prior["scale"],
                        enzyme_id=enzyme_id,
                        mic_id=mic_id,
                    )
                )
    if "transfer_constants" in parsed_toml["priors"].keys():
        for prior in parsed_toml["priors"]["transfer_constants"]:
            enzyme_id = prior["enzyme_id"]
            priors["transfer_constants"].append(
                Prior(
                    id=f"transfer_constant_{enzyme_id}",
                    location=prior["location"],
                    scale=prior["scale"],
                    enzyme_id=enzyme_id,
                )
            )
    stan_codes = get_stan_codes(kinetic_model, experiments)
    mi = MaudInput(
        kinetic_model=kinetic_model,
        priors=priors,
        stan_codes=stan_codes,
        experiments=experiments,
    )
    validation.validate_maud_input(mi)
    return mi
