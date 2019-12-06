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

from collections import defaultdict

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
    Modifier,
    Parameter,
    Prior,
    Reaction,
)
from string import ascii_lowercase
from typing import Dict

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


def load_kinetic_model_from_toml(parsed_toml: dict, model_id: str = 'mi') -> KineticModel:
    compartments = {
        c["id"]: Compartment(id=c["id"], name=c["name"], volume=c["volume"])
        for c in parsed_toml["compartments"]
    }
    metabolites = {
        m["id"] + "_" + m["compartment"]: Metabolite(
            id=m["id"],
            name=m["name"],
            balanced=m["balanced"],
            compartment=m["compartment"],
        )
        for m in parsed_toml["metabolites"]
    }
    reactions = {
        r["id"]: load_reaction_from_toml(r) for r in parsed_toml["reactions"]
    }
    return KineticModel(
        model_id=model_id,
        metabolites=metabolites,
        compartments=compartments,
        reactions=reactions
    )


def get_toml_enzyme_params(toml_enzyme: dict, stoichiometry: dict) -> Dict[str, Parameter]:
    if toml_enzyme["mechanism"] == "modular_rate_law":
        substrates = [k for k, v in stoichiometry.items() if v < 0]
        products = [k for k, v in stoichiometry.items() if v > 0]
        substrate_param_ids = [
            "K" + letter for substrate, letter in zip(substrates, ascii_lowercase)
        ]
        product_param_ids = [
            "K" + letter for product, letter in zip(products, ascii_lowercase[15:])
        ]
        param_ids = substrate_param_ids + product_param_ids + ['Kcat1']
        return {
            param_id: Parameter(param_id, toml_enzyme["id"])
            for param_id in param_ids
        }
    else:
        return {
            param_id: Parameter(param_id, toml_enzyme["id"])
            for param_id in MECHANISM_TO_PARAM_IDS[e["mechanism"]]
        }


def parse_modifiers(toml_enzyme: dict) -> tuple:
    modifier_types =  "competitive_inhibitors"
    modifiers = {}
    modifier_params = {}
    for modifier_type in ['allosteric_inhibitor', 'allosteric_activator']:
        if modifier_type + 's' in toml_enzyme.keys():
            modifier_params["transfer_constant"] = Parameter("transfer_constant", e["id"])
            for mod in toml_enzyme[modifier_type]:
                diss_const_type = "t" if modifier_type == 'allosteric_inhibitor' else "r"
                diss_const_id = "dissociation_constant_" + diss_const_type + "_" + mod
                modifiers[f"{mod}_{modifier_type}"] = Modifier(mod, modifier_type)
                modifier_params[diss_const_id] = Parameter(diss_const_id, toml_enzyme["id"], mod)
    if 'competitive_inhibitors' in toml_enzyme.keys():
        if e["mechanism"] != "modular_rate_law":
                raise ValueError(
                    """competitive inhibitors are currently
                    only supported for the mechanism 'modular_rate_law'"""
                )
        for mod in toml_enzyme["competitive_inhibitors"]:
            modifiers[f"{mod}_competitive_inhibitors"] = Modifier(mod, "competitive_inhibitor")
            inhibition_constant_id = f"inhibition_constant_{mod}"
            modifier_params[inhibition_constant_id] = Parameter(inhibition_constant_id, toml_enzyme["id"], mod)
    return modifiers, modifier_params


def load_reaction_from_toml(toml_reaction: dict) -> Reaction:
    enzymes = {}
    reversible = toml_reaction["reversible"] if "reversible" in toml_reaction.keys() else None
    is_exchange = toml_reaction["is_exchange"] if "is_exchange" in toml_reaction.keys() else None
    for e in toml_reaction["enzymes"]:
        non_modifier_params = get_toml_enzyme_params(e, toml_reaction['stoichiometry'])
        modifiers, modifier_params = parse_modifiers(e)
        params = {**non_modifier_params, **modifier_params}
        modifier_params = defaultdict()
        enzymes[e["id"]] = Enzyme(
            id=e["id"],
            name=e["name"],
            reaction_id=toml_reaction["id"],
            mechanism=e["mechanism"],
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
        for target_type in ["metabolite", "reaction"]:
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
    for exp_id, umps in parsed_toml["priors"]["unbalanced_metabolites"].items():
        for ump in umps:
            prior_id = f"{exp_id}_{ump['target_id']}"
            priors[prior_id] = Prior(
                id=prior_id,
                target_id=ump["target_id"],
                location=ump["location"],
                scale=ump["scale"],
                target_type="unbalanced_metabolite",
                experiment_id=exp_id,
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
    for exp_id, eps in parsed_toml["priors"]["enzymes"].items():
        for ep in eps:
            prior_id = f"{exp_id}_{ep['target_id']}"
            priors[prior_id] = Prior(
                id=prior_id,
                target_id=ep["target_id"],
                location=ep["location"],
                scale=ep["scale"],
                target_type="enzyme",
                experiment_id=exp_id,
            )

    mi = MaudInput(kinetic_model=kinetic_model, priors=priors, experiments=experiments)
    validation.validate_maud_input(mi)
    return mi
