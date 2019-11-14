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


def load_maud_input_from_toml(filepath: str, id: str = "mi") -> MaudInput:
    """
    Load an MaudInput object from a suitable toml file.

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
            if e["mechanism"] == "modular_rate_law":
                params = {
                    param_id["target_id"]: Parameter(param_id["target_id"], e["id"])
                    for param_id in parsed_toml["priors"]["kinetic_parameters"][e["id"]]
                    if param_id["target_id"]
                    not in [
                        "dissociation_constant_t",
                        "dissociation_constant_r",
                        "inhibition_constant",
                        "transfer_constant",
                    ]
                }
            else:
                params = {
                    param_id: Parameter(param_id, e["id"])
                    for param_id in MECHANISM_TO_PARAM_IDS[e["mechanism"]]
                }
            params["delta_g"] = Parameter("delta_g", e["id"], is_thermodynamic=True)
            modifiers = defaultdict()
            modifier_params = defaultdict()
            if any(
                [
                    x in ["allosteric_inhibitors", "allosteric_activators"]
                    for x in e.keys()
                ]
            ):
                modifier_params = {
                    "transfer_constant": Parameter("transfer_constant", e["id"])
                }
            if "allosteric_inhibitors" in e.keys():
                for inhibitor_id in e["allosteric_inhibitors"]:
                    modifiers.update(
                        {
                            f"{inhibitor_id}_allosteric_inhibitor": Modifier(
                                inhibitor_id, "allosteric_inhibitor"
                            )
                        }
                    )
                    diss_t_const_id = f"dissociation_constant_t_{inhibitor_id}"
                    modifier_params.update(
                        {
                            diss_t_const_id: Parameter(
                                diss_t_const_id, e["id"], inhibitor_id
                            )
                        }
                    )
            if "allosteric_activators" in e.keys():
                for activator_id in e["allosteric_activators"]:
                    modifiers.update(
                        {
                            f"{activator_id}_allosteric_allosteric_activator": Modifier(
                                activator_id, "allosteric_activator"
                            )
                        }
                    )
                    diss_r_const_id = f"dissociation_constant_r_{activator_id}"
                    modifier_params.update(
                        {
                            diss_r_const_id: Parameter(
                                diss_r_const_id, e["id"], activator_id
                            )
                        }
                    )
            if "competitive_inhibitors" in e.keys():
                if e["mechanism"] != "modular_rate_law":
                    raise ValueError(
                        """competitive inhibitors are currently
                        only supported for the mechanism 'modular_rate_law'"""
                    )

                for inhibitor_id in e["competitive_inhibitors"]:
                    modifiers.update(
                        {
                            f"{inhibitor_id}_competitive_inhibitors": Modifier(
                                inhibitor_id, "competitive_inhibitor"
                            )
                        }
                    )
                    inhibition_constant_id = f"inhibition_constant_{inhibitor_id}"
                    modifier_params.update(
                        {
                            inhibition_constant_id: Parameter(
                                inhibition_constant_id, e["id"], inhibitor_id
                            )
                        }
                    )
            params.update(modifier_params)
            enz = Enzyme(
                id=e["id"],
                name=e["name"],
                reaction_id=r["id"],
                mechanism=e["mechanism"],
                parameters=params,
                modifiers=modifiers,
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
    for marginal_dg in parsed_toml["priors"]["thermodynamic_parameters"][
        "marginal_dgs"
    ]:
        prior_id = f"{marginal_dg['target_id']}_delta_g"
        priors[prior_id] = Prior(
            id=prior_id,
            target_id=marginal_dg["target_id"],
            location=marginal_dg["location"],
            scale=marginal_dg["scale"],
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
