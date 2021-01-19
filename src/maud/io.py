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
import pandas as pd
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
    StanCodeSet,
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
    experiment_set = parse_experiment_set_df(pd.read_csv(experiments_path))
    prior_set = parse_prior_set_df(pd.read_csv(priors_path))
    stan_codes = get_stan_codes(kinetic_model, experiment_set)
    mi = MaudInput(
        config=config,
        kinetic_model=kinetic_model,
        experiments=experiment_set,
        priors=prior_set,
        stan_codes=stan_codes,
    )
    validation.validate_maud_input(mi)
    return mi


def parse_toml_kinetic_model(raw: dict) -> KineticModel:
    """Turn the output of toml.load into a KineticModel object.

    :param raw: Result of running toml.load on a suitable toml file

    """
    model_id = raw["model_id"] if "model_id" in raw.keys() else "mi" 
    compartments = {
        c["id"]: Compartment(id=c["id"], name=c["name"], volume=c["volume"])
        for c in raw["compartments"]
    }
    metabolites = {
        m["id"]: Metabolite(id=m["id"], name=m["name"])
        for m in raw["metabolites"]
    }
    mics = {
        f"{m['id']}_{m['compartment']}": MetaboliteInCompartment(
            id=f"{m['id']}_{m['compartment']}",
            metabolite_id=m["id"],
            compartment_id=m["compartment"],
            balanced=m["balanced"],
        )
        for m in raw["metabolites"]
    }
    reactions = {r["id"]: parse_toml_reaction(r) for r in raw["reactions"]}
    if "drains" in raw.keys():
        drains = {d["id"]: parse_toml_drain(d) for d in raw["drains"]}
    else:
        drains = None
    if "phosphorylation" in raw.keys():
        phosphorylation = {
            d["id"]: parse_toml_phosphorylation(d)
            for d in raw["phosphorylation"]
        }
    else:
        phosphorylation = None
    return KineticModel(
        model_id=model_id,
        metabolites=metabolites,
        compartments=compartments,
        mics=mics,
        reactions=reactions,
        drains=drains,
        phosphorylation=phosphorylation,
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
    if km.phosphorylation is not None:
        phos_enz_codes = codify(km.phosphorylation.keys())
    else:
        phos_enz_codes = {}
    return StanCodeSet(
        metabolite_codes=codify(km.metabolites.keys()),
        mic_codes=mic_codes,
        balanced_mic_codes=balanced_mic_codes,
        unbalanced_mic_codes=unbalanced_mic_codes,
        reaction_codes=codify(km.reactions.keys()),
        experiment_codes=codify([e.id for e in experiments.experiments]),
        enzyme_codes=codify(enzyme_ids),
        phos_enz_codes=phos_enz_codes,
        drain_codes=drain_codes,
    )


def parse_toml_drain(raw: dict) -> Drain:
    """Turn a dictionary representing a drain into a Drain object.

    :param raw: Dictionary representing a drain, typically one of
    the values of the 'drains' field in the output of toml.load.

    """

    return Drain(
        id=raw["id"],
        name=raw["name"],
        stoichiometry=raw["stoichiometry"],
    )


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
    enzymes = {}
    subunits = 1
    reversible = (
        raw["reversible"] if "reversible" in raw.keys() else None
    )
    is_exchange = (
        raw["is_exchange"] if "is_exchange" in raw.keys() else None
    )
    water_stoichiometry = (
        raw["water_stoichiometry"]
        if "water_stoichiometry" in raw.keys()
        else 0
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

        enzymes[e["id"]] = Enzyme(
            id=e["id"],
            name=e["name"],
            reaction_id=raw["id"],
            modifiers=modifiers,
            subunits=subunits,
        )
    return Reaction(
        id=raw["id"],
        name=raw["name"],
        reversible=reversible,
        is_exchange=is_exchange,
        stoichiometry=raw["stoichiometry"],
        enzymes=enzymes,
        water_stoichiometry=water_stoichiometry,
    )


def parse_experiment_set_df(raw: pd.DataFrame) -> ExperimentSet:
    """get an ExperimentSet object from a dataframe.

    :param raw: result of running pd.read_csv on a suitable file
    """
    experiment_list = []
    for id, exp_df in raw.groupby("experiment_id"):
        measurements = {"mic": {}, "flux": {}, "enzyme":{}}
        knockouts = []
        for i, row in exp_df.iterrows():
            measurement_type = row["measurement_type"]
            target_id = row["target_id"]
            if measurement_type == "knockout":
                knockouts.append(target_id)
            else:
                measurements[measurement_type][target_id] = Measurement(
                    target_id=target_id,
                    value=row["measurement"],
                    uncertainty=row["error_scale"],
                    target_type=measurement_type
                )
        experiment_list.append(
            Experiment(id=id, measurements=measurements, knockouts=knockouts)
        )
    return ExperimentSet(experiment_list)


def parse_prior_set_df(raw: pd.DataFrame) -> PriorSet:
    """get a PriorSet object from a dataframe.

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
    negative_param_types = ["formation_energy",]
    for _, row in raw.iterrows():
        id = "_".join(
            row.loc[lambda s: [isinstance(x, str) for x in s]].values
        )
        parameter_type = row["parameter_type"]
        prior_dict = row.dropna().drop("parameter_type").to_dict()
        is_non_negative = parameter_type not in negative_param_types
        priors[parameter_type + "_priors"].append(Prior(
            id=id, is_non_negative=is_non_negative, **prior_dict
        ))
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


