# Copyright (C) 2021 Novo Nordisk Foundation Center for Biosustainability,
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

"""Export yaml file combatible with yaml2sbml."""

import os

import numpy as np
import pandas as pd
from jinja2 import Template

from maud import io
from maud.sampling import get_full_stoichiometry


Template_T_met = Template(
    """{%- for met in met_array -%} ({{met[0]}} /{{met[1]}})^{{met[2]}} \
    {%- if not loop.last %} * {% endif %} {%- endfor -%}"""
)

Template_Haldane = Template(
    """{%- for Km in Km_array -%} ({{Km[0]}})^{{Km[1]}} \
    {%- if not loop.last %} * {% endif %} {%- endfor -%} / {{Keq}}"""
)

Template_Tr = Template("""{{enz}} * {{kcat}} * ({{Trf}} - {{Trr}} * {{Hal}})""")

Template_Dr = Template(
    """{%- for met in met_array -%} (1 + {{met[0]}}/{{met[1]}})^{{met[2]}} \
    {%- if not loop.last %} + {% endif %} {%- endfor -%}"""
)

Template_Drreg = Template(
    """{%- for met in met_array -%} ({{met[0]}}/{{met[1]}}) \
    {%- if not loop.last %} + {% endif %} {%- endfor -%}"""
)

Template_Allo = Template(
    """1/(1 + {{L0}}*(({{Dr}} + {{Drreg}} - 1)*\
    {{Allo_Inh}}/{{Allo_Act}})^{{Subunits}})"""
)

Template_Allo_Act_Inh = Template(
    """{%- for met in met_array -%} (1 + {{met[0]}}/{{met[1]}}) \
    {%- if not loop.last %} * {% endif %} {%- endfor -%}"""
)

Template_flux = Template("""({{Tr}})/({{Dr}} + {{Drreg}} - 1)*{{Allo}}""")

Template_drain = Template("""{{Drain}}""")


Template_yaml = Template(
    """time:
    variable: t

parameters:
    - parameterId: Zero
      nominalValue: 0
{%- for par in parameters %}
    - parameterId: {{par[0]}}
      nominalValue: {{par[1]}}
{% endfor -%}

odes:
{%- for ode in odes %}
    - stateId: {{ode[0]}}
      rightHandSide: {{ode[1]}}
      initialValue: {{ode[2]}}
    {% endfor -%}

"""
)

HERE = os.path.dirname(os.path.abspath(__file__))
relative_path_Maud_Input = "../tests/data/ecoli_small/"
relative_path_output = "output.yaml"

mi = io.load_maud_input_from_toml(os.path.join(HERE, relative_path_Maud_Input))
selected_experiment = None
if selected_experiment is None:
    selected_experiment = list(mi.stan_codes.experiment_codes.keys())[0]

kinetic_model = mi.kinetic_model
enzyme_codes = mi.stan_codes.enzyme_codes
mic_codes = mi.stan_codes.mic_codes
balanced_mic_codes = mi.stan_codes.balanced_mic_codes
met_codes = mi.stan_codes.metabolite_codes
drain_codes = mi.stan_codes.drain_codes
reaction_codes = mi.stan_codes.reaction_codes
_, _, S, _, _ = get_full_stoichiometry(
    kinetic_model, enzyme_codes, mic_codes, reaction_codes, drain_codes
)
prior_dfs = {
    prior_type: pd.DataFrame(map(lambda p: p.__dict__, priors))
    if len(priors) > 0
    else pd.DataFrame([], columns=["location", "scale", "mic_code"])
    for prior_type, priors in mi.priors.__dict__.items()
}
formation_energies = prior_dfs["formation_energy_priors"] = (
    prior_dfs["formation_energy_priors"]
    .assign(stan_code=lambda df: df["metabolite_id"].map(met_codes))
    .sort_values("stan_code")
)

flux_vector = []

for rxn in mi.kinetic_model.reactions.values():

    tmp_stoic = pd.DataFrame(rxn.stoichiometry, index=[0]).transpose()
    tmp_stoic = tmp_stoic.reset_index()
    tmp_stoic.columns = ["mic_id", "stoichiometry"]
    tmp_formation_energies = {
        mic: formation_energies[formation_energies["metabolite_id"] == mic[:-2]][
            "location"
        ].values[0]
        for mic in tmp_stoic["mic_id"].values
    }
    if rxn.water_stoichiometry:
        tmp_stoic.loc[-1] = ["h2o", rxn.water_stoichiometry]
        tmp_formation_energies.update({"h2o": -157.6})

    tmp_dG = 0
    for _, row in tmp_stoic.iterrows():
        tmp_dG += (
            0.008314
            * 298.15
            * row["stoichiometry"]
            * tmp_formation_energies[row["mic_id"]]
        )

    tmp_Keq = np.exp(tmp_dG)

    for enz_id, enz in rxn.enzymes.items():
        tmp_kms = prior_dfs["km_priors"]
        tmp_kms = tmp_kms[tmp_kms["enzyme_id"] == enz_id]
        tmp_kcat = prior_dfs["kcat_priors"]
        tmp_kcat = tmp_kcat[tmp_kcat["enzyme_id"] == enz_id]
        tmp_inh = prior_dfs["inhibition_constant_priors"]
        if len(tmp_inh) > 0:
            tmp_inh = tmp_inh[tmp_inh["enzyme_id"] == enz_id]
        tmp_allo_act = prior_dfs["relaxed_dissociation_constant_priors"]
        if len(tmp_allo_act) > 0:
            tmp_allo_act = tmp_allo_act[tmp_allo_act["enzyme_id"] == enz_id]
        tmp_allo_inh = prior_dfs["tense_dissociation_constant_priors"]
        if len(tmp_allo_inh) > 0:
            tmp_allo_inh = tmp_allo_inh[tmp_allo_inh["enzyme_id"] == enz_id]
        tmp_L = prior_dfs["transfer_constant_priors"]
        if len(tmp_L) > 0:
            tmp_L = tmp_L[tmp_L["enzyme_id"] == enz_id]

        tmp_enz_entry = tmp_kms.merge(
            tmp_stoic, left_on=["mic_id"], right_on=["mic_id"]
        )

        substrate_information = tmp_enz_entry[tmp_enz_entry["stoichiometry"] < 0]
        product_information = tmp_enz_entry[tmp_enz_entry["stoichiometry"] > 0]

        substrate_entry = zip(
            substrate_information["mic_id"],
            substrate_information["id"],
            substrate_information["stoichiometry"].abs(),
        )
        metabolite_denomintor = zip(
            tmp_enz_entry["mic_id"],
            tmp_enz_entry["id"],
            tmp_enz_entry["stoichiometry"].abs(),
        )
        product_entry = zip(
            product_information["mic_id"],
            product_information["id"],
            product_information["stoichiometry"].abs(),
        )
        haldane_entry = zip(tmp_enz_entry["id"], tmp_enz_entry["stoichiometry"])

        competitive_entry = []
        allosteric_inhibitors = []
        allosteric_activators = []

        for mod in enz.modifiers["competitive_inhibitor"]:
            competitive_entry.append(
                [mod.mic_id, tmp_inh[tmp_inh["mic_id"] == mod.mic_id].id.values[0]]
            )
        for mod in enz.modifiers["allosteric_activator"]:
            allosteric_activators.append(
                [
                    mod.mic_id,
                    tmp_allo_act[tmp_allo_act["mic_id"] == mod.mic_id].id.values[0],
                ]
            )
        for mod in enz.modifiers["allosteric_inhibitor"]:
            allosteric_inhibitors.append(
                [
                    mod.mic_id,
                    tmp_allo_inh[tmp_allo_inh["mic_id"] == mod.mic_id].id.values[0],
                ]
            )

        Drreg = Template_Drreg.render(met_array=competitive_entry)
        Trf = Template_T_met.render(met_array=substrate_entry)
        Trr = Template_T_met.render(met_array=product_entry)
        Hal = Template_Haldane.render(Km_array=haldane_entry, Keq=tmp_Keq)
        Tr = Template_Tr.render(
            enz=enz_id, kcat=tmp_kcat["id"].values[0], Trf=Trf, Trr=Trr, Hal=Hal
        )
        Dr = Template_Dr.render(met_array=metabolite_denomintor)
        if competitive_entry == []:
            Drreg = "0"
        Allo_Act = Template_Allo_Act_Inh.render(met_array=allosteric_activators)
        Allo_Inh = Template_Allo_Act_Inh.render(met_array=allosteric_inhibitors)
        if allosteric_activators == []:
            Allo_Act = "1"
        if allosteric_inhibitors == []:
            Allo_Inh = "1"
        if len(tmp_L) > 0:
            Allo = Template_Allo.render(
                L0=tmp_L["id"].values[0],
                Dr=Dr,
                Drreg=Drreg,
                Allo_Inh=Allo_Inh,
                Allo_Act=Allo_Act,
                Subunits=enz.subunits,
            )
        else:
            Allo = "1"
        flux = Template_flux.render(Tr=Tr, Dr=Dr, Drreg=Drreg, Allo=Allo)

        flux_vector.append(flux)

if any(drain_codes):
    for drain_id in kinetic_model.drains.keys():
        flux_vector.append(drain_id)

system_odes = []

for met_ix, met_vec in enumerate(S.T.values):
    if mi.kinetic_model.mics[list(mic_codes.keys())[met_ix]].balanced == 1:
        tmp_met_ode = ""
        first = 0
        for flux_ix, stoic in enumerate(met_vec):
            if stoic != 0:
                if first == 0:
                    first += 1
                    tmp_met_ode += f"({stoic}*{flux_vector[flux_ix]})"
                else:
                    tmp_met_ode += f"+({stoic}*{flux_vector[flux_ix]})"
        system_odes.append(tmp_met_ode)

par_input = []
par_input = par_input + [[par.id, par.location] for par in mi.priors.kcat_priors]
par_input = par_input + [[par.id, par.location] for par in mi.priors.km_priors]
par_input = par_input + [
    [par.id, par.location] for par in mi.priors.inhibition_constant_priors
]
par_input = par_input + [
    [par.id, par.location] for par in mi.priors.tense_dissociation_constant_priors
]
par_input = par_input + [
    [par.id, par.location] for par in mi.priors.relaxed_dissociation_constant_priors
]
par_input = par_input + [
    [par.id, par.location] for par in mi.priors.transfer_constant_priors
]
par_input = par_input + [
    [par.drain_id, par.location]
    for par in mi.priors.drain_priors
    if par.experiment_id == selected_experiment
]
par_input = par_input + [
    [par.enzyme_id, par.location]
    for par in mi.priors.enzyme_concentration_priors
    if par.experiment_id == selected_experiment
]
par_input = par_input + [
    [par.mic_id, par.location]
    for par in mi.priors.unbalanced_metabolite_priors
    if par.experiment_id == selected_experiment
]

met_init = pd.DataFrame(0.01, index=balanced_mic_codes, columns=["init"])

for exp in mi.experiments.experiments:
    if exp.id == selected_experiment:
        for meas in exp.measurements["mic"].values():
            met_init.loc[meas.target_id, "init"] = meas.value
ode_input = [
    [met, system_odes[ix], met_init.loc[met, "init"]]
    for ix, met in enumerate(balanced_mic_codes.keys())
]
yaml_input = Template_yaml.render(parameters=par_input, odes=ode_input)
with open(relative_path_output, "w") as file:
    file.writelines(yaml_input)
