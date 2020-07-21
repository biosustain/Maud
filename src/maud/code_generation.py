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

"""
Functions for generating Stan programs from MaudInput objects.

The only function that should be used outside this module is `create_stan_program`.
"""

import os
from typing import Dict, List

from jinja2 import Environment, PackageLoader, Template

from maud.data_model import MaudInput


TEMPLATE_FILES = [
    "inference_model_lower_blocks.stan",
    "functions_block.stan",
    "ode_function.stan",
    "fluxes_function.stan",
    "modular_rate_law.stan",
]


def create_stan_program(mi: MaudInput, model_type: str, time_step=0.05) -> str:
    """
    Create a stan program from an MaudInput.

    :param mi: an MaudInput object
    :param model_type: String describing the model, e.g. 'inference',
        'simulation'. So far only 'inference' is implemented.
    :param time_step: How far ahead the ode simulates
    """
    templates = get_templates()
    kinetic_model = mi.kinetic_model
    mic_codes = mi.stan_codes["metabolite_in_compartment"]
    par_codes = mi.stan_codes["kinetic_parameter"]
    balanced_codes = [
        mic_codes[mic_id] for mic_id, mic in kinetic_model.mics.items() if mic.balanced
    ]
    unbalanced_codes = [
        mic_codes[mic_id]
        for mic_id, mic in kinetic_model.mics.items()
        if not mic.balanced
    ]
    keq_position = [par_codes[par_id] for par_id in par_codes.keys() if "Keq" in par_id]
    fluxes_function = create_fluxes_function(mi, templates["fluxes_function"])
    ode_function = create_ode_function(mi, templates["ode_function"])
    functions_block = templates["functions_block"].render(
        fluxes_function=fluxes_function, ode_function=ode_function
    )
    if model_type == "inference":
        lower_blocks = templates["inference_model_lower_blocks"].render(
            balanced_codes=balanced_codes,
            unbalanced_codes=unbalanced_codes,
            Keq_position=keq_position,
        )
    else:
        raise ValueError("Model types other than 'inference' are not yet supported.")
    return functions_block + "\n" + lower_blocks


def get_templates(template_files: List[str] = TEMPLATE_FILES) -> Dict[str, Template]:
    """Load jinja templates from files.

    :param template_files: A list of paths to template files
    """

    out = {}
    env = Environment(loader=PackageLoader("maud", "stan_code"))
    for template_file in template_files:
        template_name = os.path.splitext(template_file)[0]
        out[template_name] = env.get_template(template_file)
    return out


def create_ode_function(mi: MaudInput, template: Template) -> str:
    """Get a Stan function specifying metabolite rates of change.

    :param mi: a MaudInput object
    :param template: a jinja template
    """

    kinetic_model = mi.kinetic_model
    rxn_codes = mi.stan_codes["reaction"]
    rxns = kinetic_model.reactions
    metabolite_lines = []
    for mic_id, mic in kinetic_model.mics.items():
        if mic.balanced:
            line = ""
            for rxn_id, rxn in rxns.items():
                if mic_id in rxn.stoichiometry.keys():
                    stoich = rxn.stoichiometry[mic_id]
                    prefix = "+" if stoich > 0 and line != "" else ""
                    rxn_stan = rxn_codes[rxn_id]
                    line += prefix + str(stoich) + "*fluxes[{}]".format(rxn_stan)
            metabolite_lines.append(line)
    return template.render(ode_stoic=metabolite_lines, N_flux=len(rxns))


def create_fluxes_function(mi: MaudInput, template: Template) -> str:
    """Get a Stan function that maps metabolite/parameter configurations to fluxes.

    :param mi: A MaudInput object
    :param template: A jinja template
    """
    kinetic_model = mi.kinetic_model
    modular_template = get_templates()["modular_rate_law"]
    unbalanced = [m for m in kinetic_model.mics.values() if not m.balanced]
    balanced = [m for m in kinetic_model.mics.values() if m.balanced]
    kp_codes = mi.stan_codes["kinetic_parameter"]
    enz_codes = keq_codes = mi.stan_codes["enzyme"]
    unb_ode_code = {x.id: f"p[{i+1}]" for i, x in enumerate(unbalanced)}
    bal_ode_code = {x.id: f"m[{i+1}]" for i, x in enumerate(balanced)}
    mic_ode_string = {**unb_ode_code, **bal_ode_code}

    # Add increments to codes so they can be referred to in context of stacked
    # theta vector in the Stan model. This is necessary because only one vector
    # of random variables can be passed into the algebra solver.
    enz_codes_in_theta, keq_codes_in_theta, kp_codes_in_theta = (
        {k: v + increment for k, v in original_codes.items()}
        for increment, original_codes in [
            (len(unbalanced), enz_codes),
            (len(unbalanced) + len(enz_codes), keq_codes),
            (len(unbalanced) + len(enz_codes) + len(keq_codes), kp_codes),
        ]
    )
    modular_lines = []
    free_enzyme_ratio_lines = []
    flux_lines = []
    for _, rxn in kinetic_model.reactions.items():
        irr_mass_action_rxn = 0
        for enz_id, enz in rxn.enzymes.items():
            if "irr_mass_action" == enz.mechanism:
                irr_mass_action_rxn = 1

        if irr_mass_action_rxn is 0:
            substrate_ids = [mic_id for mic_id, s in rxn.stoichiometry.items() if s < 0]
            product_ids = [mic_id for mic_id, s in rxn.stoichiometry.items() if s > 0]
            enzyme_flux_strings = []
            for enz_id, enz in rxn.enzymes.items():
                substrate_stoichiometries = [
                    [substrate_id, rxn.stoichiometry[substrate_id]]
                    for substrate_id in substrate_ids
                ]
                product_stoichiometries = [
                    [product_id, rxn.stoichiometry[product_id]]
                    for product_id in product_ids
                ]
                enz_code = enz_codes_in_theta[enz.id]
                competitor_ids = [
                    mod.mic
                    for mod in enz.modifiers.values()
                    if "competitive_inhibitor" in mod.modifier_type
                ]
                mod_substrates, mod_products, mod_competitors = get_modular_rate_codes(
                    enz_id,
                    competitor_ids,
                    substrate_stoichiometries,
                    product_stoichiometries,
                    kp_codes_in_theta,
                    mic_ode_string,
                )
                modular_line = modular_template.render(
                    enz_id=enz_id,
                    enz=enz_code,
                    Kcat1=kp_codes_in_theta[enz_id + "_" + "Kcat1"],
                    Keq=keq_codes_in_theta[enz_id],
                    substrate_list=mod_substrates,
                    product_list=mod_products,
                    competitive_inhibitor_list=mod_competitors,
                )
                modular_lines.append(modular_line)
                # make catalytic effect string
                enz_code = enz_codes[enz.id]
                catalytic_string = (
                    f"modular_rate_law(Tr_{enz_id}, Dr_{enz_id}, Dr_reg_{enz_id})"
                )
                if any([mod.allosteric for mod in enz.modifiers.values()]):
                    # make free enzyme ratio line
                    free_enzyme_ratio_line = Template(
                        "real free_enzyme_ratio_{{enzyme}} = "
                        "get_free_enzyme_ratio_{{catalytic_string}};"
                    ).render(enzyme=enz_id, catalytic_string=catalytic_string)
                    free_enzyme_ratio_lines.append(free_enzyme_ratio_line)
                    # make regulatory effect string
                    allosteric_inhibitors, allosteric_activators = (
                        {
                            mod.mic: mod
                            for mod in enz.modifiers.values()
                            if mod.modifier_type == modifier_type
                        }
                        for modifier_type in [
                            "allosteric_inhibitor",
                            "allosteric_activator",
                        ]
                    )
                    allosteric_inhibitor_codes = {
                        mod_id: mic_ode_string[mod_id]
                        for mod_id in allosteric_inhibitors.keys()
                    }
                    allosteric_activator_codes = {
                        mod_id: mic_ode_string[mod_id]
                        for mod_id in allosteric_activators.keys()
                    }
                    regulatory_string = get_regulatory_string(
                        inhibitor_codes=allosteric_inhibitor_codes,
                        activator_codes=allosteric_activator_codes,
                        num_subunits=enz.subunits,
                        param_codes=kp_codes_in_theta,
                        enzyme_name=enz.id,
                    )
                    enzyme_flux_string = catalytic_string + "*" + regulatory_string
                else:
                    enzyme_flux_string = catalytic_string
                enzyme_flux_strings.append(enzyme_flux_string)
        else:
            enz_code = enz_codes_in_theta[f"{rxn.id}"]
            cat_code = kp_codes_in_theta[f"{rxn.id}_V1"]
            enzyme_flux_strings = [f"irr_mass_action({bal_ode_code[[*rxn.stoichiometry.keys()][0]]},p[{enz_code}])"]
        flux_line = "+".join(enzyme_flux_strings)
        flux_lines.append(flux_line)
    return template.render(
        modular_coefficients=modular_lines,
        free_enzyme_ratio=free_enzyme_ratio_lines,
        fluxes=flux_lines,
    )


def get_regulatory_string(
    inhibitor_codes: Dict[str, int],
    activator_codes: Dict[str, int],
    num_subunits: int,
    param_codes: Dict[str, int],
    enzyme_name: str,
) -> str:
    """Create a call to the Stan function get_regulatory_effect.

    :param inhibitor_codes: dictionary mapping ids of inhibitor metabolites to
    integer indexes
    :param param_codes: dictionary mapping ids of parameters to integer indexes
    :param enzyme_name: string id of enzyme
    """

    regulatory_template = Template(
        """get_regulatory_effect(
        {{activator_str}},
        {{inhibitor_str}},
        free_enzyme_ratio_{{enzyme_name}},
        {{diss_r_str}},
        {{diss_t_str}},
        p[{{transfer_constant_code}}],
        {{num_subunits}}
    )""".replace(
            " ", ""
        ).replace(
            "\n", ""
        )
    )
    regulator_template = Template("{ {{-regulator_strs|join(',')-}} }")
    diss_const_template = Template("{ {{-diss_const_strs|join(',')-}} }")
    diss_t_param_codes = [
        param_codes[enzyme_name + "_dissociation_constant_t_" + inhibitor_name]
        for inhibitor_name in inhibitor_codes.keys()
    ]
    diss_r_param_codes = [
        param_codes[enzyme_name + "_dissociation_constant_r_" + activator_name]
        for activator_name in activator_codes.keys()
    ]
    transfer_constant_code = param_codes[enzyme_name + "_transfer_constant"]
    inhibitor_strs = [c for c in inhibitor_codes.values()]
    activator_strs = [c for c in activator_codes.values()]
    diss_t_strs = [f"p[{c}]" for c in diss_t_param_codes]
    diss_r_strs = [f"p[{c}]" for c in diss_r_param_codes]
    diss_t_str = diss_const_template.render(diss_const_strs=diss_t_strs)
    diss_r_str = diss_const_template.render(diss_const_strs=diss_r_strs)
    inhibitor_str = regulator_template.render(regulator_strs=inhibitor_strs)
    activator_str = regulator_template.render(regulator_strs=activator_strs)

    if not activator_codes:
        activator_str = "empty_array"
        diss_r_str = "empty_array"

    if not inhibitor_codes:
        inhibitor_str = "empty_array"
        diss_t_str = "empty_array"

    return regulatory_template.render(
        activator_str=activator_str,
        inhibitor_str=inhibitor_str,
        enzyme_name=enzyme_name,
        diss_r_str=diss_r_str,
        diss_t_str=diss_t_str,
        transfer_constant_code=transfer_constant_code,
        num_subunits=num_subunits,
    )


def get_modular_rate_codes(
    enz_id: str,
    competitor_ids: List[List],
    substrate_info: List[List],
    product_info: List[List],
    par_codes: Dict[str, str],
    mic_ode_string: Dict[str, str],
) -> List[List[int]]:
    """Get codes that can be put into the modular rate law jinja template.

    The function returns a list containing lists substrate_input and
    product_input. Each of these is a list containing lists with the form
    [mic_code, param_code, stoic] for each substrate and product.

    :param rxn_id: id of the reaction
    :param substrate_info: list containing lists with the form [mic_id, stoic]
    where mic_id is a string and stoic is a float
    :param product_info: list containing lists with the form [mic_id, stoic]
    where mic_id is a string and stoic is a float
    :param par_codes: dictionary mapping parameter ids to integers
    :param mic_ode_string: dictionary mapping metabolite-in-compartment ids to
    stan ode code strings
    """

    substrate_keys = ["a", "b", "c", "d"]
    product_keys = ["p", "q", "r", "s"]
    substrate_input = []
    product_input = []
    competitor_input = []
    for info, keys in zip(
        [substrate_info, product_info], [substrate_keys, product_keys]
    ):
        for i, (mic_id, stoic) in enumerate(info):
            param_id = enz_id + "_K" + keys[i]
            param_code = par_codes[param_id]
            if stoic < 0:
                substrate_input.append([mic_ode_string[mic_id], param_code, stoic])
            elif stoic > 0:
                product_input.append([mic_ode_string[mic_id], param_code, stoic])
    for comp in competitor_ids:
        competitor_code = mic_ode_string[comp]
        param_id = enz_id + "_inhibition_constant_" + comp
        competitor_parameter = par_codes[param_id]
        competitor_input.append([competitor_code, competitor_parameter])
    return [substrate_input, product_input, competitor_input]
