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

from string import ascii_lowercase
from typing import Dict, List, Tuple

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
    Parameter,
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
    "irr_mass_action": ["V1"]
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
    kinetic_parameter_ids = [
        f"{enzyme_id}_{parameter_id}"
        for reaction in km.reactions.values()
        for enzyme_id, enzyme in reaction.enzymes.items()
        for parameter_id in enzyme.parameters.keys()
    ]
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
        "kinetic_parameter": codify(kinetic_parameter_ids),
        "reaction": codify(km.reactions.keys()),
        "experiment": codify(experiments.keys()),
        "enzyme": codify(enzyme_ids),
    }


def get_non_modifier_params(
    enzyme_id: str, enzyme_mechanism: str, stoichiometry: dict
) -> Dict[str, Parameter]:
    """Get the non-modifer parameters for an enzyme.

    :param enzyme_id: String identifying the enzyme
    :param enzyme_mechanism: String identifying the enzyme's mechanism,
    e.g. 'ordered_bibi'
    :param stoichiometry: Dictionary representing the reaction's stoichiometry. The
    keys are metabolite ids and the values are stoichiometries.

    """
    if enzyme_mechanism == "modular_rate_law":
        substrates = [k for k, v in stoichiometry.items() if v < 0]
        products = [k for k, v in stoichiometry.items() if v > 0]
        substrate_param_ids = [
            "K" + letter for substrate, letter in zip(substrates, ascii_lowercase)
        ]
        product_param_ids = [
            "K" + letter for product, letter in zip(products, ascii_lowercase[15:])
        ]
        param_ids = substrate_param_ids + product_param_ids + ["Kcat1"]
        return {param_id: Parameter(param_id, enzyme_id) for param_id in param_ids}
    else:
        return {
            param_id: Parameter(param_id, enzyme_id)
            for param_id in MECHANISM_TO_PARAM_IDS[enzyme_mechanism]
        }


def parse_modifiers(
    modifier_type: str, modifiers: List[str], enzyme_id: str, param_prefix: str
) -> Tuple:
    """Turn unformatted modifier info into dictionaries of Modifier and Parameter objects.

    :param modifier_type: String identifying the type of modifier. One of
    'allosteric_inhibitor', 'allosteric_activator', 'competitive_inhibitor'

    :param modifier: List of strings identifying modifying metabolites.

    :param enzyme id: String identifying the enzyme that is modified.

    :param_prefix: String that goes before the metabolite id in order to
    identify the metabolite-specific parameters.

    """
    modifier_dict = {}
    modifier_params = {}
    for modifier in modifiers:
        modifier_dict[f"{modifier}_{modifier_type}"] = Modifier(modifier, modifier_type)
        param_id = f"{param_prefix}_{modifier}"
        modifier_params[param_id] = Parameter(param_id, enzyme_id, modifier)
    if modifier_type in ["allosteric_inhibitor", "allosteric_activator"]:
        modifier_params["transfer_constant"] = Parameter("transfer_constant", enzyme_id)
    return modifier_dict, modifier_params


def load_reaction_from_toml(toml_reaction: dict) -> Reaction:
    """Turn a dictionary representing a reaction into a Reaction object.

    :param toml_reaction: Dictionary representing a reaction, typically one of
    the values of the 'reactions' field in the output of toml.load.

    """
    enzymes = {}
    reversible = (
        toml_reaction["reversible"] if "reversible" in toml_reaction.keys() else None
    )
    is_exchange = (
        toml_reaction["is_exchange"] if "is_exchange" in toml_reaction.keys() else None
    )
    for e in toml_reaction["enzymes"]:
        non_modifier_params = get_non_modifier_params(
            e["id"], e["mechanism"], toml_reaction["stoichiometry"]
        )
        modifiers = {}
        modifier_params = {}
        if "competitive_inhibitors" in e.keys() & e["mechanism"] != "modular_rate_law":
            raise ValueError(
                """competitive inhibitors are currently
                only supported for the mechanism 'modular_rate_law'"""
            )
        for modifier_type, param_prefix in [
            ("allosteric_activator", "dissociation_constant_r"),
            ("allosteric_inhibitor", "dissociation_constant_t"),
            ("competitive_inhibitor", "inhibition_constant"),
        ]:
            if modifier_type + "s" in e.keys():
                modifiers_of_this_type, modifier_params_of_this_type = parse_modifiers(
                    modifier_type, e[modifier_type + "s"], e["id"], param_prefix
                )
                modifiers.update(modifiers_of_this_type)
                modifier_params.update(modifier_params_of_this_type)
        params = {**non_modifier_params, **modifier_params}
        if 'subunits' in e.keys():
            subunits = e['subunits']
        else:
            subunits = None
        enzymes[e["id"]] = Enzyme(
            id=e["id"],
            name=e["name"],
            reaction_id=toml_reaction["id"],
            mechanism=e["mechanism"],
            subunits=subunits,
            parameters=params,
            modifiers=modifiers,
        )
    return Reaction(
        id=toml_reaction["id"],
        name=toml_reaction["name"],
        reversible=reversible,
        is_exchange=is_exchange,
        stoichiometry=toml_reaction["stoichiometry"],
        enzymes=enzymes,
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
    priors = {}
    for formation_energy_prior in parsed_toml["priors"]["thermodynamic_parameters"][
        "formation_energies"
    ]:
        prior_id = f"{formation_energy_prior['target_id']}_formation_energy"
        priors[prior_id] = Prior(
            id=prior_id,
            target_id=formation_energy_prior["target_id"],
            location=formation_energy_prior["location"],
            scale=formation_energy_prior["scale"],
            target_type="thermodynamic_parameter",
        )
    for enz_id, kpps in parsed_toml["priors"]["kinetic_parameters"].items():
        for kpp in kpps:
            prior_id = f"{enz_id}_{kpp['target_id']}"
            if "metabolite" in kpp.keys():
                prior_id += "_" + kpp["metabolite"]
            priors[prior_id] = Prior(
                id=prior_id,
                target_id=kpp["target_id"],
                location=kpp["location"],
                scale=kpp["scale"],
                target_type="kinetic_parameter",
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
