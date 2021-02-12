from maud import io, sampling
from jinja2 import Environment, PackageLoader, Template
import pandas as pd
import os

HERE = os.path.dirname(os.path.abspath(__file__))

Template_T_met = Template("""{%- for met in met_array -%} ({{met[0]}}/{{met[1]}})^{{met[2]}} {%- if not loop.last %} * {% endif %} {%- endfor -%}""")

Template_Haldane =  Template("""{%- for Km in Km_array -%} ({{Km[0]}})^{{Km[1]}} {%- if not loop.last %} * {% endif %} {%- endfor -%} / {{Keq}}""")

Template_Tr = Template("""{{enz}} * {{kcat}} * ({{Trf}} - {{Trr}} * {{Hal}})""")

Template_Dr = Template("""{%- for met in met_array -%} (1 + {{met[0]}}/{{met[1]}})^{{met[2]}} {%- if not loop.last %} + {% endif %} {%- endfor -%}""")

Template_Drreg = Template("""{%- for met in met_array -%} ({{met[0]}}/{{met[1]}}) {%- if not loop.last %} + {% endif %} {%- endfor -%}""")

Template_Allo = Template("""1/(1 + {{L0}}*(({{Dr}} + {{Drreg}} - 1)*{{Allo_Inh}}/{{Allo_Act}})^{{Subunits}})""")

Template_Allo_Act_Inh = Template("""{%- for met in met_array -%} (1 + {{met[0]}}/{{met[1]}}) {%- if not loop.last %} * {% endif %} {%- endfor -%}""")

Template_flux = Template("""({{Tr}})/({{Dr}} + {{Drreg}} - 1)*{{Allo}}""")



Template_yaml = Template("""time:
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
      initialValue: 1
    {% endfor -%}

""")

mi = io.load_maud_input_from_toml(os.path.join(HERE, "../../tests/data/linear.toml"))
kinetic_model = mi.kinetic_model
enzyme_codes = mi.stan_codes['enzyme']
mic_codes = mi.stan_codes["metabolite_in_compartment"]

flux_vector = []

for rxn in mi.kinetic_model.reactions.values():

    tmp_stoich = pd.DataFrame(rxn.stoichiometry, index=[0]).transpose()
    tmp_stoich = tmp_stoich.reset_index()
    tmp_stoich.columns = ["mic_id", "stoichiometry"]

    tmp_Keq = 1

    for enz_id, enz in rxn.enzymes.items():
        tmp_kms = {
            prior.id: {
                "mic_id": prior.mic_id,
                "prior_mean": prior.location,
            }
            for prior_id, prior_entry in mi.priors.items() if prior_id == "kms"
            for prior in prior_entry if prior.enzyme_id == enz_id
        }

        tmp_kcat = {
            prior.id: {
                "prior_mean": prior.location,
            }
            for prior_id, prior_entry in mi.priors.items() if prior_id == "kcats"
            for prior in prior_entry if prior.enzyme_id == enz_id
        }

        tmp_kms_df = pd.DataFrame.from_dict(tmp_kms).transpose().reset_index()
        tmp_kms_df.columns = ["parameter_id", "mic_id", "prior_mean"]

        tmp_enz_entry = tmp_kms_df.merge(tmp_stoich, left_on=['mic_id'], right_on=['mic_id'])

        substrate_information = tmp_enz_entry[tmp_enz_entry['stoichiometry'] < 0]
        product_information = tmp_enz_entry[tmp_enz_entry['stoichiometry'] > 0]

        substrate_entry = zip(substrate_information['mic_id'], substrate_information['parameter_id'], substrate_information['stoichiometry'].abs())
        metabolite_denomintor = zip(tmp_enz_entry['mic_id'], tmp_enz_entry['parameter_id'], tmp_enz_entry['stoichiometry'].abs())
        product_entry = zip(product_information['mic_id'], product_information['parameter_id'], product_information['stoichiometry'].abs())
        haldane_entry = zip(tmp_enz_entry['parameter_id'], tmp_enz_entry['stoichiometry'])
        

        
        Trf = Template_T_met.render(met_array=substrate_entry)
        Trr = Template_T_met.render(met_array=product_entry)
        Hal = Template_Haldane.render(Km_array=haldane_entry, Keq=1)
        Tr  = Template_Tr.render(enz=enz_id,
                                 kcat=list(tmp_kcat.keys())[0],
                                 Trf=Trf,
                                 Trr=Trr,
                                 Hal=Hal)
        Dr = Template_Dr.render(met_array=metabolite_denomintor)
        flux = Template_flux.render(Tr=Tr,
                                    Dr=Dr,
                                    Drreg="0")

        flux_vector.append(flux)

S = sampling.get_full_stoichiometry(kinetic_model,
                                    enzyme_codes,
                                    mic_codes)
system_odes = []

for met_ix, met_vec in enumerate(S.T.values):
    tmp_met_ode = ''
    first = 0
    if mi.kinetic_model.mics[list(mic_codes.keys())[met_ix]].balanced == 1:
        for flux_ix, stoic in enumerate(met_vec):
            if stoic != 0:
                if first == 0:
                    first += 1
                    tmp_met_ode += f'({stoic}*{flux_vector[flux_ix]})'
                else:
                    tmp_met_ode += f'+({stoic}*{flux_vector[flux_ix]})' 
    else:
        tmp_met_ode = "Zero"

    system_odes.append(tmp_met_ode)

ode_input = [[met, system_odes[ix]] for ix, met in enumerate(mic_codes.keys())]
par_input = [[par.id, par.location] for prior in mi.priors.values()
                                    for par in prior]
par_input.append(['r1', 1])
par_input.append(['r2', 1])
par_input.append(['r3', 1])

yaml_input = Template_yaml.render(parameters=par_input,
                                  odes=ode_input)
with open("test.yaml", 'w') as file:
    file.writelines(yaml_input)



