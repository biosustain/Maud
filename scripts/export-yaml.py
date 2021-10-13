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

import argparse
import os

import numpy as np
from jinja2 import Template

from maud.analysis import load_infd
from maud.cli import get_inits_from_draw
from maud.io import load_maud_input_from_toml
from maud.sampling import get_stoichiometry


HELP_MSG = """
This script is to export a system of ODEs described using a yaml format.
The output can be parsed to yaml2sbml to convert the parameter set
to an sbml which can be simulated using software such as COPASI.
You are required to have a simulation/sample of your model already in an
output, this is done so that you can specify a simulation you want to inspect
further.

If you are interested in a specific experiment please define this as well,
otherwise it will take the first one defined, as it's not set up to do
multiple experiments simultaneously.

For further information on the workflow please refer to:
`https://github.com/yaml2sbml-dev/yaml2sbml`.

Prior to converting to sbml you will need to use the cmd:
`pip install yaml2sbml`
"""

HERE = os.path.dirname(os.path.abspath(__file__))

Template_T_met = Template(
    """{%- for met in met_array -%} ({{met[0]}} /{{met[1]}})^({{met[2]}}) \
    {%- if not loop.last %} * {% endif %} {%- endfor -%}"""
)

Template_Haldane = Template(
    """{%- for Km in Km_array -%} ({{Km[0]}})^{{Km[1]}} \
    {%- if not loop.last %} * {% endif %} {%- endfor -%} / {{Keq}}"""
)

Template_Tr = Template("""{{enz}} * {{kcat}} * ({{Trf}} - {{Trr}} * {{Hal}})""")

Template_Tr_irr = Template("""{{enz}} * {{kcat}} * {{Trf}}""")

Template_Dr = Template(
    """{%- for met in sub_array -%} (1 + {{met[0]}}/{{met[1]}})^({{met[2]}}) \
    {%- if not loop.last %} * {% endif %} {%- endfor -%} + \
    {%- for met in prod_array -%} (1 + {{met[0]}}/{{met[1]}})^({{met[2]}}) \
    {%- if not loop.last %} * {% endif %} {%- endfor -%}"""
)

Template_Dr_irr = Template(
    """{%- for met in sub_array -%} (1 + {{met[0]}}/{{met[1]}})^({{met[2]}}) \
    {%- if not loop.last %} * {% endif %} {%- endfor -%}"""
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

Template_drain = Template(
    """{{drain}}*{%- for met in sub_array -%} ({{met}} /({{met}} + 0.000001)) \
    {%- if not loop.last %} * {% endif %} {%- endfor -%}"""
)


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


def get_csvs(data_path):
    """Return list of csv files in maud_output samples folder."""
    return [
        os.path.join(data_path, "samples", f)
        for f in os.listdir(os.path.join(data_path, "samples"))
        if f.endswith(".csv")
    ]


def main():
    """Run the script."""
    parser = argparse.ArgumentParser(description=HELP_MSG)
    parser.add_argument(
        "maud_output_dir", type=str, nargs=1, help="A path to Maud output directory"
    )
    parser.add_argument(
        "--chain", default=0, help="Chain number to export parameter values for"
    )
    parser.add_argument(
        "--draw", default=0, help="Draw number to export parameter values for"
    )
    parser.add_argument("--warmup", default=0, help="If draw is in warmup phase or not")
    parser.add_argument(
        "--yaml_output", default="output.yaml", help="Output of constructed yaml"
    )
    parser.add_argument(
        "--selected_experiment", default=None, help="Experiment parameters exported"
    )
    args = parser.parse_args()
    maud_output_dir = args.maud_output_dir[0]
    chain = args.chain
    draw = args.draw
    warmup = args.warmup
    yaml_output = os.path.join(HERE, args.yaml_output)
    csvs = get_csvs(maud_output_dir)
    mi = load_maud_input_from_toml(os.path.join(maud_output_dir, "user_input"))
    infd = load_infd(csvs, mi)
    selected_experiment = None
    if selected_experiment is None:
        selected_experiment = list(mi.stan_coords.experiments)[0]

    # Defining Stoichiometry
    S = get_stoichiometry(mi)

    # Selecting a set of parameters from a previous run
    par_values = get_inits_from_draw(infd, mi, chain, draw, warmup)
    par_input = []

    # Selecting measurements
    conc_measurements = mi.measurements.yconc
    balanced_conc_values = conc_measurements.loc[selected_experiment]
    balanced_mic_values = {
        mic.id: balanced_conc_values.loc[mic.id]["measurement"]
        if mic.id in balanced_conc_values.index
        else 0.001
        for mic in mi.kinetic_model.mics
        if mic.balanced
    }

    # Experiment specific parameters
    exp_values = par_values[par_values["experiment_id"] == selected_experiment]
    conc_values = exp_values[exp_values["parameter_name"] == "conc_unbalanced"]
    enz_values = exp_values[exp_values["parameter_name"] == "conc_enzyme"]
    drain_values = exp_values[exp_values["parameter_name"] == "drain"]
    for mic in mi.kinetic_model.mics:
        if mic.balanced is False:
            par_input.append(
                [
                    f"m{mic.id}",
                    list(conc_values[conc_values["mic_id"] == mic.id]["value"])[0],
                ]
            )
    for rxn in mi.kinetic_model.reactions:
        for enz in rxn.enzymes:
            par_input.append(
                [
                    f"e{enz.id}",
                    list(enz_values[enz_values["enzyme_id"] == enz.id]["value"])[0],
                ]
            )
    for rxn in mi.kinetic_model.reactions:
        if rxn.reaction_mechanism == "drain":
            par_input.append(
                [
                    f"r{rxn.id}",
                    list(drain_values[drain_values["drain_id"] == rxn.id]["value"])[0],
                ]
            )

    # Metabolite gibbs energies
    dgfs = par_values[par_values["parameter_name"] == "dgf"]

    flux_dict = {}

    for rxn in mi.kinetic_model.reactions:
        # calculating the gibbs energy of reaction
        # accounting for water
        if rxn.reaction_mechanism == "reversible_modular_rate_law":
            tmp_dg = 0
            for mic_id, stoic in rxn.stoichiometry.items():
                met_id = next(
                    filter(lambda mic: mic.id == mic_id, mi.kinetic_model.mics)
                ).metabolite_id
                met_dgf = list(dgfs[dgfs["metabolite_id"] == met_id]["value"])[0]
                tmp_dg += stoic * met_dgf
            if rxn.water_stoichiometry:
                tmp_dg += rxn.water_stoichiometry * -157.6
            tmp_Keq = np.exp(tmp_dg / (-0.008314 * 298.15))

        for enz in rxn.enzymes:
            tmp_enz_pars = par_values[par_values["enzyme_id"] == enz.id]
            tmp_kms = tmp_enz_pars.loc[tmp_enz_pars["parameter_name"] == "km"]
            par_input += [
                [f"km_{row['enzyme_id']}_{row['mic_id']}", row["value"]]
                for _, row in tmp_kms.iterrows()
            ]
            tmp_kcats = tmp_enz_pars.loc[tmp_enz_pars["parameter_name"] == "kcat"]
            par_input += [
                [f"kcat_{row['enzyme_id']}", row["value"]]
                for _, row in tmp_kcats.iterrows()
            ]
            tmp_kis = tmp_enz_pars.loc[tmp_enz_pars["parameter_name"] == "ki"]
            par_input += [
                [f"ki_{row['enzyme_id']}_{row['mic_id']}", row["value"]]
                for _, row in tmp_kis.iterrows()
            ]
            tmp_aas = tmp_enz_pars.loc[tmp_enz_pars["parameter_name"] == "diss_r"]
            par_input += [
                [f"aa_{row['enzyme_id']}_{row['mic_id']}", row["value"]]
                for _, row in tmp_aas.iterrows()
            ]
            tmp_ais = tmp_enz_pars.loc[tmp_enz_pars["parameter_name"] == "diss_t"]
            par_input += [
                [f"ai_{row['enzyme_id']}_{row['mic_id']}", row["value"]]
                for _, row in tmp_ais.iterrows()
            ]
            tmp_transfer_constants = tmp_enz_pars.loc[
                tmp_enz_pars["parameter_name"] == "transfer_constant"
            ]
            par_input += [
                [f"transfer_constant_{row['enzyme_id']}", row["value"]]
                for _, row in tmp_transfer_constants.iterrows()
            ]

            substrate_list = [
                f"m{mic}" for mic, stoic in rxn.stoichiometry.items() if stoic < 0
            ]
            product_list = [
                f"m{mic}" for mic, stoic in rxn.stoichiometry.items() if stoic > 0
            ]
            mic_list = [f"m{mic}" for mic, _ in rxn.stoichiometry.items()]

            substrate_entry = list(
                zip(
                    substrate_list,
                    [f"km_{enz.id}_{mic[1:]}" for mic in substrate_list],
                    [np.abs(rxn.stoichiometry[mic[1:]]) for mic in substrate_list],
                )
            )

            product_entry = list(
                zip(
                    product_list,
                    [f"km_{enz.id}_{mic[1:]}" for mic in product_list],
                    [np.abs(rxn.stoichiometry[mic[1:]]) for mic in product_list],
                )
            )

            haldane_entry = list(
                zip(
                    [f"km_{enz.id}_{mic[1:]}" for mic in mic_list],
                    [rxn.stoichiometry[mic[1:]] for mic in mic_list],
                )
            )

            competitive_entry = []
            allosteric_inhibitors = []
            allosteric_activators = []

            for mod in enz.modifiers["competitive_inhibitor"]:
                competitive_entry.append(
                    [f"m{mod.mic_id}", f"ki_{enz.id}_{mod.mic_id}"]
                )
            for mod in enz.modifiers["allosteric_activator"]:
                allosteric_activators.append(
                    [f"m{mod.mic_id}", f"aa_{enz.id}_{mod.mic_id}"]
                )
            for mod in enz.modifiers["allosteric_inhibitor"]:
                allosteric_inhibitors.append(
                    [f"m{mod.mic_id}", f"ai_{enz.id}_{mod.mic_id}"]
                )

            if rxn.reaction_mechanism == "reversible_modular_rate_law":
                Trf = Template_T_met.render(met_array=substrate_entry)
                Trr = Template_T_met.render(met_array=product_entry)
                Hal = Template_Haldane.render(Km_array=haldane_entry, Keq=tmp_Keq)
                Tr = Template_Tr.render(
                    enz=f"e{enz.id}", kcat=f"kcat_{enz.id}", Trf=Trf, Trr=Trr, Hal=Hal
                )
                Dr = Template_Dr.render(
                    sub_array=substrate_entry, prod_array=product_entry
                )

            elif rxn.reaction_mechanism == "irreversible_modular_rate_law":
                Trf = Template_T_met.render(met_array=substrate_entry)
                Tr = Template_Tr_irr.render(
                    enz=f"e{enz.id}", kcat=f"kcat_{enz.id}", Trf=Trf
                )
                Dr = Template_Dr_irr.render(sub_array=substrate_entry)

            Drreg = Template_Drreg.render(sub_array=substrate_entry)
            if competitive_entry == []:
                Drreg = "0"
            Allo_Act = Template_Allo_Act_Inh.render(met_array=allosteric_activators)
            Allo_Inh = Template_Allo_Act_Inh.render(met_array=allosteric_inhibitors)
            if allosteric_activators == []:
                Allo_Act = "1"
            if allosteric_inhibitors == []:
                Allo_Inh = "1"
            if any([allosteric_inhibitors, allosteric_activators]):
                Allo = Template_Allo.render(
                    L0=f"transfer_constant_{enz.id}",
                    Dr=Dr,
                    Drreg=Drreg,
                    Allo_Inh=Allo_Inh,
                    Allo_Act=Allo_Act,
                    Subunits=enz.subunits,
                )
            else:
                Allo = "1"
            flux = Template_flux.render(Tr=Tr, Dr=Dr, Drreg=Drreg, Allo=Allo)
            flux_dict[enz.id] = flux
        if rxn.reaction_mechanism == "drain":
            substrate_list = [
                f"m{mic}" for mic, stoic in rxn.stoichiometry.items() if stoic < 0
            ]
            if substrate_list == []:
                substrate_list = [1]
            flux = Template_drain.render(drain=f"r{rxn.id}", sub_array=substrate_list)
            flux_dict[rxn.id] = flux

    system_odes = {}
    for mic in mi.kinetic_model.mics:
        if mic.balanced is True:
            tmp_met_ode = ""
            first = 0
            for edge in mi.stan_coords.edges:
                if S.loc[mic.id, edge] != 0:
                    if first == 0:
                        first += 1
                        tmp_met_ode += f"({S.loc[mic.id, edge]}*{flux_dict[edge]})"
                    else:
                        tmp_met_ode += f"+({S.loc[mic.id, edge]}*{flux_dict[edge]})"
            system_odes[mic.id] = tmp_met_ode
    ode_input = [
        [f"m{mic.id}", system_odes[mic.id], balanced_mic_values[mic.id]]
        for mic in mi.kinetic_model.mics
        if mic.balanced is True
    ]
    yaml_input = Template_yaml.render(parameters=par_input, odes=ode_input)
    with open(yaml_output, "w") as file:
        file.writelines(yaml_input)


if __name__ == "__main__":
    main()
