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

"""Functions for validating Maud objects."""

from maud import data_model


def validate_maud_input(mi: data_model.MaudInput):
    """Check that priors, experiments and kinetic model are consistent."""
    km = mi.kinetic_model
    km_unb_mets = [mic.id for mic in km.mics.values() if not mic.balanced]
    km_balanced_mets = [mic.id for mic in km.mics.values() if mic.balanced]
    km_mets_no_compartments = list(set([met.id for met in km.metabolites.values()]))
    km_rxns = [rxn.id for rxn in km.reactions.values()]
    km_pars = [
        p
        for rxn in km.reactions.values()
        for enz in rxn.enzymes.values()
        for p in enz.parameters.values()
    ]
    km_kps = [p.enzyme_id + "_" + p.id for p in km_pars if not p.is_thermodynamic]
    prior_kps, prior_tps = (
        [pid for pid, p in mi.priors.items() if p.target_type == t]
        for t in ["kinetic_parameter", "thermodynamic_parameter"]
    )
    prior_enzs, prior_unb_mets = (
        [
            (p.experiment_id, p.target_id)
            for p in mi.priors.values()
            if p.target_type == t
        ]
        for t in ["enzyme", "unbalanced_metabolite"]
    )
    for kp in prior_kps:
        if kp not in km_kps:
            raise ValueError(
                f"kinetic parameter {kp} has a prior but is not in the kinetic model."
            )
    for kp in km_kps:
        if kp not in prior_kps:
            raise ValueError(
                f"kinetic parameter {kp} is in the kinetic model but has no prior."
            )
    for tp in prior_tps:
        tp_met = tp.replace("_formation_energy", "")
        if tp_met not in km_mets_no_compartments:
            raise ValueError(
                f"Metabolite {tp_met} has a formation energy prior but is not in the "
                "kinetic model."
            )
    for km_met in km_mets_no_compartments:
        if km_met + "_formation_energy" not in prior_tps:
            raise ValueError(
                f"Metabolite {km_met} is in the kinetic model but has no formation"
                " energy prior."
            )
    for exp in mi.experiments.values():
        for meas in exp.measurements["metabolite"].values():
            if meas.target_id not in km_unb_mets + km_balanced_mets:
                raise ValueError(
                    f"metabolite {meas.target_id} is measured in experiment {exp.id}"
                    "but is not in the kinetic model {mi.kinetic_model.model_id}."
                )
        for meas in exp.measurements["reaction"].values():
            if meas.target_id not in km_rxns:
                raise ValueError(
                    f"reaction {meas.target_id} is measured in experiment {exp.id}"
                    "but is not in the kinetic model {mi.kinetic_model.model_id}."
                )
