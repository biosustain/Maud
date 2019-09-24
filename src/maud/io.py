# Copyright (c) 2019, Novo Nordisk Foundation Center for Biosustainability, Technical University of Denmark.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from collections import defaultdict

import toml

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


MECHANISM_TO_PARAM_IDS = {"uniuni": ["Keq", "Kcat1", "Kcat2", "Ka"]}


def load_maud_input_from_toml(filepath: str, id: str = "mi") -> MaudInput:
    """
    Load an MaudInput object from a suitable toml file
    
    :param filepath: path to a toml file
    :param id: id for the output object

    """
    kinetic_model = KineticModel(id)
    parsed_toml = toml.load(filepath)
    for c in parsed_toml["compartments"]:
        cmp = Compartment(id=c["id"], name=c["name"], volume=c["volume"])
        kinetic_model.compartments.update({cmp.id: cmp})
    for m in parsed_toml["metabolites"]:
        met = Metabolite(
            id=m["id"],
            name=m["name"],
            balanced=m["balanced"],
            compartment=m["compartment"],
        )
        kinetic_model.metabolites.update({met.id: met})
    for r in parsed_toml["reactions"]:
        rxn_enzymes = {}
        for e in r["enzymes"]:
            params = {
                param_id: Parameter(param_id, e["id"])
                for param_id in MECHANISM_TO_PARAM_IDS[e["mechanism"]]
            }
            allosteric_inhibitors = defaultdict()
            if "allosteric_inhibitors" in e.keys():
                allosteric_params = {
                    "transfer_constant": Parameter("transfer_constant", e["id"])
                }
                for inhibitor_id in e["allosteric_inhibitors"]:
                    allosteric_inhibitors.update(
                        {inhibitor_id: Modifier(inhibitor_id, "allosteric_inhibitor")}
                    )
                    diss_t_const_id = f"dissociation_constant_t_{inhibitor_id}"
                    allosteric_params.update(
                        {
                            diss_t_const_id: Parameter(
                                diss_t_const_id, e["id"], inhibitor_id
                            )
                        }
                    )
                params.update(allosteric_params)
            enz = Enzyme(
                id=e["id"],
                name=e["name"],
                reaction_id=r["id"],
                mechanism=e["mechanism"],
                parameters=params,
                modifiers=allosteric_inhibitors,
            )
            rxn_enzymes.update({enz.id: enz})
        rxn = Reaction(
            id=r["id"],
            name=r["name"],
            reversible=r["reversible"] if "reversible" in r.keys() else None,
            is_exchange=r["is_exchange"] if "is_exchange" in r.keys() else None,
            stoichiometry=r["stoichiometry"],
            enzymes=rxn_enzymes,
        )
        kinetic_model.reactions.update({rxn.id: rxn})
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

    return MaudInput(
        kinetic_model=kinetic_model, priors=priors, experiments=experiments
    )
