import os
import pandas as pd
from enzymekat._data_model import (
    EnzymeKatInput,
    KineticModel,
    Metabolite,
    Reaction,
    Modifier,
    Parameter,
    Prior,
    Experiment,
    Measurement,
    Compartment
)
from jinja2 import Template, Environment, PackageLoader, select_autoescape
from textwrap import dedent
from typing import Dict, List, Iterable

JINJA_TEMPLATES = [
    'inference_model_template.stan',
    'simulation_model_template.stan',
    'functions_block.stan',
    'ode_function.stan',
    'fluxes_function.stan',
    'steady_state_function.stan'
]
MECHANISM_TEMPLATES = {
    'uniuni': Template(
        "uniuni(m[{{S0}}], m[{{P0}}], e[{{enz}}]*p[{{Kcat1}}], e[{{enz}}]*p[{{Kcat2}}], p[{{Ka}}], p[{{Keq}}])"
    )
}


def codify(l: Iterable[str]):
    return dict(zip(l, range(1, len(l) + 1)))


def get_templates(template_files=JINJA_TEMPLATES):
    out = {}
    env = Environment(
        loader=PackageLoader('enzymekat', 'stan_code'),
        autoescape=select_autoescape(['stan'])
    )
    for template_file in template_files:
        template_name = os.path.splitext(template_file)[0]
        out[template_name] = env.get_template(template_file)
    return out


def create_stan_program(eki: EnzymeKatInput, model_type: str) -> str:
    """
    Create a stan program from an EnzymeKatInput.

    :param eki: an EnzymeKatInput object
    :param model_type: one of 'inference_model', 'simulation_model'
    """
    templates = get_templates()
    kinetic_model = eki.kinetic_model
    fluxes_function = create_fluxes_function(
        kinetic_model, templates['fluxes_function']
    )
    ode_function = create_ode_function(
        kinetic_model, templates['ode_function']
    )
    steady_state_function = create_steady_state_function(
        kinetic_model, templates['steady_state_function']
    )
    functions_block = templates['functions_block'].render(
        fluxes_function=fluxes_function,
        ode_function=ode_function,
        steady_state_function=steady_state_function
    )
    lower_blocks = templates[model_type].render()
    return functions_block + '\n' + lower_blocks


def create_ode_function(kinetic_model: KineticModel, template: Template) -> str:
    rxns = kinetic_model.reactions
    rxn_id_to_stan = codify(rxns.keys())
    metabolite_lines = []
    for met_id, met in kinetic_model.metabolites.items():
        if not met.balanced:
            line = '0'
        else:
            line = ''
            for rxn_id, rxn in rxns.items():
                if met_id in rxn.stoichiometry.keys():
                    stoich = rxn.stoichiometry[met_id]
                    prefix = '+' if stoich > 0 and line != '' else ''
                    rxn_stan = rxn_id_to_stan[rxn_id]
                    line += prefix + str(stoich) + '*fluxes[{}]'.format(rxn_stan)
        metabolite_lines.append(line)
    return template.render(ode_stoic=metabolite_lines, N_flux=len(rxns))
                    

def create_steady_state_function(
        kinetic_model: KineticModel,
        template: Template,
        time_step: float = 0.05
) -> str:
    mets = kinetic_model.metabolites
    met_codes = dict(zip(mets.keys(), range(1, len(mets) + 1)))
    balanced_codes = [
        met_codes[met_id] for met_id, met in mets.items() if met.balanced
    ]
    unbalanced_codes = [
        met_codes[met_id] for met_id, met in mets.items() if not met.balanced
    ]
    return template.render(
        N_balanced=len(balanced_codes),
        N_unbalanced=len(unbalanced_codes),
        time_step=time_step,
        balanced_codes=balanced_codes,
        unbalanced_codes=unbalanced_codes
    )


def create_Kip_ordered_unibi_line(param_codes: dict, rxn_id: str) -> str:
    template = Template(
        "real {{rxn_id}}_Kip = get_Kip_ordered_unibi({{Keq}}, {{Ka}}, {{Kp}}, {{Kcat1}}, {{Kcat2}});"
    )
    return template.render(
        rxn_id=rxn_id,
        Keq=param_codes[rxn_id + '_Keq'],
        Kia=param_codes[rxn_id + 'Kia'],
        Kq=param_codes[rxn_id + 'Kq'],
        Kcat1=param_codes[rxn_id + 'Kcat1'],
        Kcat2=param_codes[rxn_id + 'Kcat2'],
    )


def create_Kiq_ordered_unibi_line(param_codes: dict, rxn_id: str) -> str:
    template = Template(
        "real {{rxn_id}}_Kiq = get_Kiq_ordered_unibi({{Keq}}, {{Ka}}, {{Kp}}, {{Kcat1}}, {{Kcat2}});"
    )
    return template.render(
        rxn_id=rxn_id,
        Keq=param_codes[rxn_id + 'Keq'],
        Kia=param_codes[rxn_id + 'Kia'],
        Kq=param_codes[rxn_id + 'Kq'],
        Kcat1=param_codes[rxn_id + 'Kcat1'],
        Kcat2=param_codes[rxn_id + 'Kcat2'],
    )


def create_fluxes_function(kinetic_model: KineticModel, template: Template) -> str:
    mechanism_to_haldane_functions = {
        'ordered_unibi': [
            create_Kip_ordered_unibi_line,
            create_Kiq_ordered_unibi_line
        ],
    }
    par_codes = get_parameter_codes(kinetic_model)
    met_codes = get_metabolite_codes(kinetic_model)
    enz_codes = get_enzyme_codes(kinetic_model)
    haldane_lines = []
    free_enzyme_ratio_lines = []
    flux_lines = []
    for rxn_id, rxn in kinetic_model.reactions.items():
        substrate_ids = [met_id for met_id, s in rxn.stoichiometry.items() if s > 0]
        substrate_codes = {
            'S' + str(i): met_codes[met_id] for i, met_id in enumerate(substrate_ids)
        }
        product_ids = [met_id for met_id, s in rxn.stoichiometry.items() if s < 0]
        product_codes = {
            'P' + str(i): met_codes[met_id] for i, met_id in enumerate(product_ids)
        }
        enzyme_flux_strings = []
        for enz_id, enz in rxn.enzymes.items():
            # make haldane relationship lines if necessary
            if enz.mechanism in mechanism_to_haldane_functions.keys():
                haldane_functions = mechanism_to_haldane_functions[rxn.mechanism]
                for haldane_function in haldane_functions:
                    haldane_line = haldane_function(par_codes, enz.id)
                    haldane_lines.append(haldane_line)
            # make catalytic effect string
            enz_param_codes = {k[len(enz_id)+1:]: v for k, v in par_codes.items() if enz_id in k}
            enz_code = enz_codes[enz.id]
            mechanism_args = {**substrate_codes, **product_codes, **{'enz': enz_code}, **enz_param_codes}
            catalytic_string = MECHANISM_TEMPLATES[enz.mechanism].render(mechanism_args)
            if any(enz.modifiers):
                # make free enzyme ratio line
                free_enzyme_ratio_line = Template(
                    "real free_enzyme_ratio_{{enzyme}} = get_free_enzyme_ratio_{{catalytic_string}};"
                ).render(enzyme=enz_id, catalytic_string=catalytic_string)
                free_enzyme_ratio_lines.append(free_enzyme_ratio_line)
                # make regulatory effect string
                allosteric_inhibitors = {
                    mod_id: mod for mod_id, mod in enz.modifiers.items()
                    if mod.modifier_type == 'allosteric_inhibitor'
                }
                allosteric_inhibitor_codes = {
                    mod_id: met_codes[mod_id] for mod_id in allosteric_inhibitors.keys()
                }
                regulatory_string = get_regulatory_string(
                    allosteric_inhibitor_codes, par_codes, enz.id
                )
                enzyme_flux_string = catalytic_string + '*' + regulatory_string
            else:
                enzyme_flux_string = catalytic_string
            enzyme_flux_strings.append(enzyme_flux_string)
        flux_line = '+'.join(enzyme_flux_strings)
        flux_lines.append(flux_line)
    return template.render(
        haldanes=haldane_lines,
        free_enzyme_ratio=free_enzyme_ratio_lines,
        fluxes=flux_lines
    )


def get_regulatory_string(inhibitor_codes, param_codes, enzyme_name):
    regulatory_template = Template("""get_regulatory_effect(
        {},
        {{inhibitor_str}},
        free_enzyme_ratio_{{enzyme_name}},
        {},
        {{diss_t_str}},
        p[{{transfer_constant_code}}]
    )""".replace(' ', '').replace('\n', ''))
    inhibitor_template = Template("{ {{-inhibitor_strs|join(',')-}} }")
    diss_t_template = Template("{ {{-diss_t_strs|join(',')-}} }")
    diss_t_param_codes = [
        param_codes[enzyme_name + '_dissociation_constant_t_' + inhibitor_name]
        for inhibitor_name in inhibitor_codes.keys()
    ]
    transfer_constant_code = param_codes[enzyme_name + '_transfer_constant']
    inhibitor_strs = [f'm[{c}]' for c in inhibitor_codes.values()]
    diss_t_strs = [f'p[{c}]' for c in diss_t_param_codes]
    diss_t_str = diss_t_template.render(diss_t_strs=diss_t_strs)
    inhibitor_str = inhibitor_template.render(inhibitor_strs=inhibitor_strs)
    return regulatory_template.render(
        inhibitor_str=inhibitor_str,
        enzyme_name=enzyme_name,
        diss_t_str=diss_t_str,
        transfer_constant_code=transfer_constant_code
    )


def get_metabolite_codes(kinetic_model: KineticModel) -> Dict[str, int]:
    return codify(kinetic_model.metabolites.keys())


def get_enzyme_codes(kinetic_model: KineticModel) -> Dict[str, int]:
    enzyme_ids = []
    for rxn_id, rxn in kinetic_model.reactions.items():
        for enz_id, _ in rxn.enzymes.items():
            enzyme_ids.append(enz_id)
    return codify(enzyme_ids)


def get_parameter_codes(kinetic_model: KineticModel) -> Dict[str, int]:
    parameter_ids = []
    for rxn_id, rxn in kinetic_model.reactions.items():
        for enz_id, enz in rxn.enzymes.items():
            for par_id, par in enz.parameters.items():
                parameter_ids.append(enz_id + '_' + par_id)
    return codify(parameter_ids)
