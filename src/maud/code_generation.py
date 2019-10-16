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

from maud.data_model import KineticModel, MaudInput
from maud.utils import codify


TEMPLATE_FILES = [
    "inference_model_lower_blocks.stan",
    "functions_block.stan",
    "ode_function.stan",
    "fluxes_function.stan",
    "steady_state_function.stan",
    "modular_rate_law.stan",
]
MECHANISM_TEMPLATES = {
    "uniuni": Template(
        "uniuni(m[{{S0}}],m[{{P0}}],p[{{enz}}]*p[{{Kcat1}}],"
        "p[{{enz}}]*p[{{Kcat2}}],p[{{Ka}}],p[{{Keq}}])"
    ),
    "ordered_unibi": Template(
        "ordered_unibi(m[{{S0}}],m[{{P0}}],m[{{P1}}],"
        "p[{{enz}}]*p[{{Kcat1}}],p[{{enz}}]*p[{{Kcat2}}],"
        "p[{{Ka}}],p[{{Kp}}],p[{{Kq}}],"
        "p[{{Kia}}],{{enz_id}}_Kip,{{enz_id}}_Kiq],"
        "p[{{Keq}}])"
    ),
    "ordered_bibi": Template(
        "ordered_bibi(m[{{S0}}],m[{{S1}}],m[{{P0}}],m[{{P1}}],"
        "p[{{enz}}]*p[{{Kcat1}}],p[{{enz}}]*p[{{Kcat2}}],"
        "p[{{Ka}}],p[{{Kb}}],p[{{Kp}}],p[{{Kq}}],"
        "{{enz_id}}_Kia,p[{{Kib}}],{{enz_id}}_Kip,p[{{Kiq}}],"
        "p[{{Keq}}])"
    ),
    "pingpong": Template(
        "pingpong(m[{{S0}}],m[{{S1}}],m[{{P0}}],m[{{P1}}],"
        "p[{{enz}}]*p[{{Kcat1}}],p[{{enz}}]*p[{{Kcat2}}],"
        "p[{{Ka}}],p[{{Kb}}],p[{{Kp}}],p[{{Kq}}],"
        "p[{{Kia}}],p[{{Kib}}],{{enz_id}}_Kip,p[{{Kiq}}],"
        "p[{{Keq}}])"
    ),
    "ordered_terbi": Template(
        "ordered_terbi(m[{{S0}}],m[{{S1}}],m[{{S2}}],m[{{P0}}],m[{{P1}}],"
        "p[{{enz}}]*p[{{Kcat1}}],p[{{enz}}]*p[{{Kcat2}}],"
        "p[{{Ka}}],p[{{Kb}}],p[{{Kc}}],{{enz_id}}_Kp,p[{{Kq}}],"
        "p[{{Kia}}],p[{{Kib}}],p[{{Kic}}],{{enz_id}}_Kip,p[{{Kiq}}],"
        "p[{{Keq}}])"
    ),
    "modular_rate_law": Template("modular_rate_law({{Tr}},{{Dr}})"),
}


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
    met_codes = get_metabolite_codes(kinetic_model)
    par_codes = get_parameter_codes(kinetic_model)
    balanced_codes = [
        met_codes[met_id]
        for met_id, met in kinetic_model.metabolites.items()
        if met.balanced
    ]
    unbalanced_codes = [
        met_codes[met_id]
        for met_id, met in kinetic_model.metabolites.items()
        if not met.balanced
    ]
    keq_position = [par_codes[par_id] for par_id in par_codes.keys() if "Keq" in par_id]
    fluxes_function = create_fluxes_function(
        kinetic_model, templates["fluxes_function"]
    )
    ode_function = create_ode_function(kinetic_model, templates["ode_function"])
    steady_state_function = create_steady_state_function(
        kinetic_model, templates["steady_state_function"], time_step
    )
    functions_block = templates["functions_block"].render(
        fluxes_function=fluxes_function,
        ode_function=ode_function,
        steady_state_function=steady_state_function,
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


def create_ode_function(kinetic_model: KineticModel, template: Template) -> str:
    """Get a Stan function specifying metabolite rates of change.

    :param kinetic_model: a KineticModel object
    :param template: a jinja template
    """

    rxns = kinetic_model.reactions
    rxn_id_to_stan = codify(rxns.keys())
    metabolite_lines = []
    for met_id, met in kinetic_model.metabolites.items():
        if not met.balanced:
            line = "0"
        else:
            line = ""
            for rxn_id, rxn in rxns.items():
                if met_id in rxn.stoichiometry.keys():
                    stoich = rxn.stoichiometry[met_id]
                    prefix = "+" if stoich > 0 and line != "" else ""
                    rxn_stan = rxn_id_to_stan[rxn_id]
                    line += prefix + str(stoich) + "*fluxes[{}]".format(rxn_stan)
        metabolite_lines.append(line)
    return template.render(ode_stoic=metabolite_lines, N_flux=len(rxns))


def create_steady_state_function(
    kinetic_model: KineticModel, template: Template, time_step: float
) -> str:
    """Get a Stan function for finding the steady state of the system.

    This is a function that can be input to Stan's algebra solver, whose return
    value is zero when the system is at steady state.

    :param kinetic_model: A KineticModel object
    :param template: A jinja template
    :param time_step: A number specifying how far into the future the ode
    solver should simulate the system in order to find an evolved state to
    compare with the initial state
    """

    mets = kinetic_model.metabolites
    met_codes = dict(zip(mets.keys(), range(1, len(mets) + 1)))
    balanced_codes = [met_codes[met_id] for met_id, met in mets.items() if met.balanced]
    unbalanced_codes = [
        met_codes[met_id] for met_id, met in mets.items() if not met.balanced
    ]
    return template.render(
        N_balanced=len(balanced_codes),
        N_unbalanced=len(unbalanced_codes),
        time_step=time_step,
        balanced_codes=balanced_codes,
        unbalanced_codes=unbalanced_codes,
    )


def create_Kip_ordered_unibi_line(param_codes: dict, rxn_id: str) -> str:
    """Generate Kip for mechanism ordered_unibi using the Haldane relationship."""
    template = Template(
        "real {{rxn_id}}_Kip = get_Kip_ordered_unibi(p[{{Keq}}], p[{{Ka}}], p[{{Kp}}], "
        "p[{{Kcat1}}], p[{{Kcat2}}]);"
    )
    return template.render(
        rxn_id=rxn_id,
        Keq=param_codes[rxn_id + "_" + "Keq"],
        Kia=param_codes[rxn_id + "_" + "Kia"],
        Kq=param_codes[rxn_id + "_" + "Kq"],
        Kcat1=param_codes[rxn_id + "_" + "Kcat1"],
        Kcat2=param_codes[rxn_id + "_" + "Kcat2"],
    )


def create_Kiq_ordered_unibi_line(param_codes: dict, rxn_id: str) -> str:
    """Generate Kiq for mechanism ordered_unibi using the Haldane relationship."""
    template = Template(
        "real {{rxn_id}}_Kiq = get_Kiq_ordered_unibi(p[{{Keq}}], p[{{Ka}}], p[{{Kp}}], "
        "p[{{Kcat1}}], p[{{Kcat2}}]);"
    )
    return template.render(
        rxn_id=rxn_id,
        Keq=param_codes[rxn_id + "_" + "Keq"],
        Kia=param_codes[rxn_id + "_" + "Kia"],
        Kq=param_codes[rxn_id + "_" + "Kq"],
        Kcat1=param_codes[rxn_id + "_" + "Kcat1"],
        Kcat2=param_codes[rxn_id + "_" + "Kcat2"],
    )


def create_Kia_ordered_bibi_line(param_codes: dict, rxn_id: str) -> str:
    """Generate Kia for mechanism ordered_bibi using the Haldane relationship."""
    template = Template(
        "real {{rxn_id}}_Kia = get_Kia_ordered_bibi(p[{{Keq}}], p[{{Ka}}], p[{{Kq}}], "
        "p[{{Kib}}], p[{{Kcat1}}], p[{{Kcat2}}]);"
    )
    return template.render(
        rxn_id=rxn_id,
        Keq=param_codes[rxn_id + "_" + "Keq"],
        Ka=param_codes[rxn_id + "_" + "Ka"],
        Kq=param_codes[rxn_id + "_" + "Kq"],
        Kib=param_codes[rxn_id + "_" + "Kib"],
        Kcat1=param_codes[rxn_id + "_" + "Kcat1"],
        Kcat2=param_codes[rxn_id + "_" + "Kcat2"],
    )


def create_Kip_ordered_bibi_line(param_codes: dict, rxn_id: str) -> str:
    """Generate Kip for mechanism ordered_bibi using the Haldane relationship."""
    template = Template(
        "real {{rxn_id}}_Kip = get_Kip_ordered_bibi(p[{{Keq}}], p[{{Kb}}], p[{{Kp}}], "
        "p[{{Kiq}}], p[{{Kcat1}}], p[{{Kcat2}}]);"
    )
    return template.render(
        rxn_id=rxn_id,
        Keq=param_codes[rxn_id + "_" + "Keq"],
        Kb=param_codes[rxn_id + "_" + "Kb"],
        Kp=param_codes[rxn_id + "_" + "Kp"],
        Kiq=param_codes[rxn_id + "_" + "Kiq"],
        Kcat1=param_codes[rxn_id + "_" + "Kcat1"],
        Kcat2=param_codes[rxn_id + "_" + "Kcat2"],
    )


def create_Kip_pingpong_line(param_codes: dict, rxn_id: str) -> str:
    """Generate Kip for mechanism pingpong using the Haldane relationship."""
    template = Template(
        "real {{rxn_id}}_Kip = get_Kip_pingpong({{Keq}}, {{Ka}}, {{Kb}}, {{Kq}},"
        " {{Kcat1}}, {{Kcat2}});"
    )
    return template.render(
        rxn_id=rxn_id,
        Keq=param_codes[rxn_id + "_" + "Keq"],
        Ka=param_codes[rxn_id + "_" + "Ka"],
        Kb=param_codes[rxn_id + "_" + "Kb"],
        Kq=param_codes[rxn_id + "_" + "Kq"],
        Kcat1=param_codes[rxn_id + "_" + "Kcat1"],
        Kcat2=param_codes[rxn_id + "_" + "Kcat2"],
    )


def create_Kip_ordered_terbi_line(param_codes: dict, rxn_id: str) -> str:
    """Generate Kip for mechanism ordered_terbi using the Haldane relationship."""
    template = Template(
        "real {{rxn_id}}_Kip = get_Kip_ordered_terbi({{Keq}}, {{Kc}}, {{Kia}}, "
        "{{Kib}}, {{Kiq}}, {{Kcat1}}, {{Kcat2}});"
    )
    return template.render(
        rxn_id=rxn_id,
        Keq=param_codes[rxn_id + "_" + "Keq"],
        Kc=param_codes[rxn_id + "_" + "Kc"],
        Kia=param_codes[rxn_id + "_" + "Kia"],
        Kib=param_codes[rxn_id + "_" + "Kib"],
        Kiq=param_codes[rxn_id + "_" + "Kiq"],
        Kcat1=param_codes[rxn_id + "_" + "Kcat1"],
        Kcat2=param_codes[rxn_id + "_" + "Kcat2"],
    )


def create_Kp_ordered_terbi_line(param_codes: dict, rxn_id: str) -> str:
    """Generate Kp for mechanism ordered_terbi using the Haldane relationship."""
    template = Template(
        "real {{rxn_id}}_Kp = get_Kiq_ordered_terbi({{Keq}}, {{Kia}}, {{Kib}}, "
        "{{Kic}}, {{Kiq}});"
    )
    return template.render(
        rxn_id=rxn_id,
        Keq=param_codes[rxn_id + "_" + "Keq"],
        Kia=param_codes[rxn_id + "_" + "Kia"],
        Kib=param_codes[rxn_id + "_" + "Kib"],
        Kic=param_codes[rxn_id + "_" + "Kic"],
        Kiq=param_codes[rxn_id + "_" + "Kiq"],
    )


def create_fluxes_function(kinetic_model: KineticModel, template: Template) -> str:
    """Get a Stan function that maps metabolite/parameter configurations to fluxes.

    :param kinetic_model: A KineticModel object
    :param template: A jinja template
    """
    modular_template = get_templates()["modular_rate_law"]
    mechanism_to_haldane_functions = {
        "ordered_unibi": [create_Kip_ordered_unibi_line, create_Kiq_ordered_unibi_line],
        "ordered_bibi": [create_Kia_ordered_bibi_line, create_Kip_ordered_bibi_line],
        "pingpong": [create_Kip_pingpong_line],
        "ordered_terbi": [create_Kip_ordered_terbi_line, create_Kp_ordered_terbi_line],
    }
    par_codes = get_parameter_codes(kinetic_model)
    met_codes = get_metabolite_codes(kinetic_model)
    enz_codes = get_enzyme_codes(kinetic_model)
    par_codes_in_enz_context = {k: v + len(enz_codes) for k, v in par_codes.items()}
    haldane_lines = []
    modular_lines = []
    free_enzyme_ratio_lines = []
    flux_lines = []
    for _, rxn in kinetic_model.reactions.items():
        substrate_ids = [met_id for met_id, s in rxn.stoichiometry.items() if s < 0]
        substrate_codes = {
            "S" + str(i): met_codes[met_id] for i, met_id in enumerate(substrate_ids)
        }
        product_ids = [met_id for met_id, s in rxn.stoichiometry.items() if s > 0]
        product_codes = {
            "P" + str(i): met_codes[met_id] for i, met_id in enumerate(product_ids)
        }
        enzyme_flux_strings = []
        for enz_id, enz in rxn.enzymes.items():
            # make haldane relationship lines if necessary
            if enz.mechanism in mechanism_to_haldane_functions.keys():
                haldane_functions = mechanism_to_haldane_functions[enz.mechanism]
                for haldane_function in haldane_functions:
                    haldane_line = haldane_function(par_codes_in_enz_context, enz.id)
                    haldane_lines.append(haldane_line)
            substrate_stoichiometries = [
                [substrate_id, rxn.stoichiometry[substrate_id]]
                for substrate_id in substrate_ids
            ]
            product_stoichiometries = [
                [product_id, rxn.stoichiometry[product_id]]
                for product_id in product_ids
            ]
            # make modular rate law if necessary
            if enz.mechanism == "modular_rate_law":
                enz_code = enz_codes[enz.id]
                substrate_block, product_block = get_modular_rate_codes(
                    enz_id,
                    substrate_stoichiometries,
                    product_stoichiometries,
                    par_codes_in_enz_context,
                    met_codes,
                )
                modular_line = modular_template.render(
                    enz_id=enz_id,
                    enz=enz_code,
                    Kcat1=par_codes_in_enz_context[enz_id + "_" + "Kcat1"],
                    Keq=par_codes_in_enz_context[enz_id + "_" + "Keq"],
                    substrate_list=substrate_block,
                    product_list=product_block,
                )
                modular_lines.append(modular_line)
            # make catalytic effect string
            enz_param_codes = {
                k[len(enz_id) + 1 :]: v
                for k, v in par_codes_in_enz_context.items()
                if enz_id in k
            }
            enz_code = enz_codes[enz.id]
            if enz.mechanism == "modular_rate_law":
                mechanism_args = {
                    "Tr": ("Tr_{}").format(enz_id),
                    "Dr": ("Dr_{}").format(enz_id),
                }
            else:
                mechanism_args = {
                    **substrate_codes,
                    **product_codes,
                    **{"enz": enz_code},
                    **enz_param_codes,
                    "enz_id": enz_id,
                }

            catalytic_string = MECHANISM_TEMPLATES[enz.mechanism].render(mechanism_args)
            if any(enz.modifiers):
                # make free enzyme ratio line
                free_enzyme_ratio_line = Template(
                    "real free_enzyme_ratio_{{enzyme}} = "
                    "get_free_enzyme_ratio_{{catalytic_string}};"
                ).render(enzyme=enz_id, catalytic_string=catalytic_string)
                free_enzyme_ratio_lines.append(free_enzyme_ratio_line)
                # make regulatory effect string
                allosteric_inhibitors = {
                    mod_id: mod
                    for mod_id, mod in enz.modifiers.items()
                    if mod.modifier_type == "allosteric_inhibitor"
                }
                allosteric_inhibitor_codes = {
                    mod_id: met_codes[mod_id] for mod_id in allosteric_inhibitors.keys()
                }
                regulatory_string = get_regulatory_string(
                    allosteric_inhibitor_codes, par_codes_in_enz_context, enz.id
                )
                enzyme_flux_string = catalytic_string + "*" + regulatory_string
            else:
                enzyme_flux_string = catalytic_string
            enzyme_flux_strings.append(enzyme_flux_string)
        flux_line = "+".join(enzyme_flux_strings)
        flux_lines.append(flux_line)
    return template.render(
        haldanes=haldane_lines,
        modular_coefficients=modular_lines,
        free_enzyme_ratio=free_enzyme_ratio_lines,
        fluxes=flux_lines,
    )


def get_regulatory_string(
    inhibitor_codes: Dict[str, int], param_codes: Dict[str, int], enzyme_name: str
) -> str:
    """Create a call to the Stan function get_regulatory_effect.

    :param inhibitor_codes: dictionary mapping ids of inhibitor metabolites to
    integer indexes
    :param param_codes: dictionary mapping ids of parameters to integer indexes
    :param enzyme_name: string id of enzyme
    """

    regulatory_template = Template(
        """get_regulatory_effect(
        empty_array,
        {{inhibitor_str}},
        free_enzyme_ratio_{{enzyme_name}},
        empty_array,
        {{diss_t_str}},
        p[{{transfer_constant_code}}]
    )""".replace(
            " ", ""
        ).replace(
            "\n", ""
        )
    )
    inhibitor_template = Template("{ {{-inhibitor_strs|join(',')-}} }")
    diss_t_template = Template("{ {{-diss_t_strs|join(',')-}} }")
    diss_t_param_codes = [
        param_codes[enzyme_name + "_dissociation_constant_t_" + inhibitor_name]
        for inhibitor_name in inhibitor_codes.keys()
    ]
    transfer_constant_code = param_codes[enzyme_name + "_transfer_constant"]
    inhibitor_strs = [f"m[{c}]" for c in inhibitor_codes.values()]
    diss_t_strs = [f"p[{c}]" for c in diss_t_param_codes]
    diss_t_str = diss_t_template.render(diss_t_strs=diss_t_strs)
    inhibitor_str = inhibitor_template.render(inhibitor_strs=inhibitor_strs)
    return regulatory_template.render(
        inhibitor_str=inhibitor_str,
        enzyme_name=enzyme_name,
        diss_t_str=diss_t_str,
        transfer_constant_code=transfer_constant_code,
    )


def get_modular_rate_codes(
    rxn_id: str,
    substrate_info: List[List],
    product_info: List[List],
    par_codes: Dict[str, int],
    met_codes: Dict[str, int],
):
    """Get codes that can be put into the modular rate law jinja template.

    The function returns a list containing lists substrate_input and
    product_input. Each of these is a list containing lists with the form
    [met_code, param_code, stoic] for each substrate and product.

    :param rxn_id: id of the reaction
    :param substrate_info: list containing lists with the form [met_id, stoic]
    where met_id is a string and stoic is a float
    :param product_info: list containing lists with the form [met_id, stoic]
    where met_id is a string and stoic is a float
    :param par_codes: dictionary mapping parameter ids to integers
    :param met_codes: dictionary mapping metabolite ids to integers
    """

    substrate_keys = ["a", "b", "c", "d"]
    product_keys = ["p", "q", "r", "s"]
    substrate_input = []
    product_input = []
    for info, keys in zip(
        [substrate_info, product_info], [substrate_keys, product_keys]
    ):
        for i, (met_id, stoic) in enumerate(info):
            param_id = rxn_id + "_K" + keys[i]
            param_code = par_codes[param_id]
            met_code = met_codes[met_id]
            if stoic < 0:
                substrate_input.append([met_code, param_code, stoic])
            elif stoic > 0:
                product_input.append([met_code, param_code, stoic])
    return [substrate_input, product_input]


def get_metabolite_codes(kinetic_model: KineticModel) -> Dict[str, int]:
    """Get a dictionary mapping metabolite ids to integer indexes.

    :param kinetic_model: A KineticModel object
    """

    return codify(kinetic_model.metabolites.keys())


def get_enzyme_codes(kinetic_model: KineticModel) -> Dict[str, int]:
    """Get a dictionary mapping enzyme ids to integer indexes.

    :param kinetic_model: A KineticModel object
    """

    enzyme_ids = []
    for _, rxn in kinetic_model.reactions.items():
        for enz_id, _ in rxn.enzymes.items():
            enzyme_ids.append(enz_id)
    return codify(enzyme_ids)


def get_parameter_codes(kinetic_model: KineticModel) -> Dict[str, int]:
    """Get a dictionary mapping parameter ids to integer indexes.

    :param kinetic_model: A KineticModel object
    """

    parameter_ids = []
    for _, rxn in kinetic_model.reactions.items():
        for enz_id, enz in rxn.enzymes.items():
            for par_id, _ in enz.parameters.items():
                parameter_ids.append(enz_id + "_" + par_id)
    return codify(parameter_ids)
