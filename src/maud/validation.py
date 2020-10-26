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
    model = mi.kinetic_model
    model_unb_mics = [mic.id for mic in model.mics.values() if not mic.balanced]
    model_balanced_mics = [mic.id for mic in model.mics.values() if mic.balanced]
    model_metabolites = list(set([met.id for met in model.metabolites.values()]))
    model_rxns = [rxn.id for rxn in model.reactions.values()]
    model_kms = [
        f"km_{enz.id}_{mic_id}"
        for rxn in model.reactions.values()
        for enz in rxn.enzymes.values()
        for mic_id in rxn.stoichiometry.keys()
    ]
    model_kcats = [
        f"kcat_{enz.id}"
        for rxn in model.reactions.values()
        for enz in rxn.enzymes.values()
    ]
    model_formation_energies = [
        f"formation_energy_{met_id}" for met_id in model_metabolites
    ]
    model_kis = [
        f"ki_{enz.id}_{modifier.mic_id}"
        for rxn in model.reactions.values()
        for enz in rxn.enzymes.values()
        for modifier in enz.modifiers["competitive_inhibitor"]
    ]
    prior_kms = [p.id for p in mi.priors["kms"]]
    prior_kcats = [p.id for p in mi.priors["kcats"]]
    prior_kis = [p.id for p in mi.priors["inhibition_constants"]]
    prior_formation_energies = [p.id for p in mi.priors["formation_energies"]]
    prior_drains = [p.id for p in mi.priors["drains"]]
    for model_pars, prior_pars in zip(
        [model_kms, model_kcats, model_formation_energies, model_kis],
        [prior_kms, prior_kcats, prior_formation_energies, prior_kis],
    ):
        for prior_par in prior_pars:
            msg = f"{prior_par} is in the priors but not the kinetic model."
            if prior_par not in model_pars:
                raise ValueError(msg)
        for model_par in model_pars:
            msg = f"{model_par} is in the kinetic model but not the priors."
            if model_par not in prior_pars:
                raise ValueError(msg)
    for exp in mi.experiments.values():
        for meas in exp.measurements["metabolite"].values():
            if meas.target_id not in model_unb_mics + model_balanced_mics:
                raise ValueError(
                    f"metabolite {meas.target_id} is measured in experiment {exp.id}"
                    "but is not in the kinetic model {mi.kinetic_model.model_id}."
                )
        for meas in exp.measurements["reaction"].values():
            if meas.target_id not in model_rxns:
                raise ValueError(
                    f"reaction {meas.target_id} is measured in experiment {exp.id}"
                    "but is not in the kinetic model {mi.kinetic_model.model_id}."
                )
        if mi.kinetic_model.drains is not None:
            for drain_id in mi.kinetic_model.drains:
                for drain in mi.priors["drains"]:
                    if (drain_id != drain.drain_id) & (exp.id != drain.experiment_id):
                        raise ValueError(
                            f"drain {drain_id} was not included in experiment {exp.id}."
                            "A drain prior is required for every drain and each experiment."
                        )
