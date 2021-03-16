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
    model_unb_mics = [mic.id for mic in model.mics if not mic.balanced]
    model_balanced_mics = [mic.id for mic in model.mics if mic.balanced]
    model_metabolites = list(set(met.id for met in model.metabolites))
    model_rxns = [rxn.id for rxn in model.reactions]
    model_kms = [
        f"km_{mic_id}_{enz.id}"
        for rxn in model.reactions
        for enz in rxn.enzymes
        for mic_id in rxn.stoichiometry.keys()
    ]
    experiment_ids = list(set(m.experiment_id for m in mi.measurements))
    model_kcats = [f"kcat_{enz.id}" for rxn in model.reactions for enz in rxn.enzymes]
    model_formation_energies = [
        f"formation_energy_{met_id}" for met_id in model_metabolites
    ]
    model_kis = [
        f"inhibition_constant_{modifier.mic_id}_{enz.id}"
        for rxn in model.reactions
        for enz in rxn.enzymes
        for modifier in enz.modifiers["competitive_inhibitor"]
    ]
    prior_kms = [p.id for p in mi.priors.km_priors]
    prior_kcats = [p.id for p in mi.priors.kcat_priors]
    prior_kis = [p.id for p in mi.priors.inhibition_constant_priors]
    prior_formation_energies = [p.id for p in mi.priors.formation_energy_priors]
    for model_pars, prior_pars in zip(
        [model_kms, model_kcats, model_formation_energies, model_kis],
        [prior_kms, prior_kcats, prior_formation_energies, prior_kis],
    ):
        for prior_par in prior_pars:
            msg = f"{prior_par} is in the priors but not the kinetic model."
            if prior_par not in model_pars:
                print(model_pars)
                raise ValueError(msg)
        for model_par in model_pars:
            msg = f"{model_par} is in the kinetic model but not the priors."
            if model_par not in prior_pars:
                raise ValueError(msg)
    for meas in mi.measurements:
        if meas.target_type == "mic":
            if meas.target_id not in model_unb_mics + model_balanced_mics:
                raise ValueError(
                    f"metabolite {meas.target_id} is measured in experiment"
                    "{meas.experiment_id} but is not in the kinetic model "
                    "{mi.kinetic_model.model_id}."
                )
        elif meas.target_type == "flux":
            if meas.target_id not in model_rxns:
                raise ValueError(
                    f"reaction {meas.target_id} is measured in experiment"
                    "{meas.experiment_id} but is not in the kinetic model"
                    "{mi.kinetic_model.model_id}."
                )
    for drain in mi.kinetic_model.drains:
        for experiment_id in experiment_ids:
            if not any(
                p.drain_id == drain.id and p.experiment_id == experiment_id
                for p in mi.priors.drain_priors
            ):
                raise ValueError(
                    f"Drain {drain.id} has no prior for experiment {experiment_id}."
                )
