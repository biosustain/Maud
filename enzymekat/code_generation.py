import os
import pandas as pd
from enzymekat.data_model import EnzymeKatData


TEMPLATE_RELATIVE_PATHS = {
    'inference': 'stan_code/inference_model_template.stan',
    'simulation': 'stan_code/simulation_model_template.stan',
}


def create_stan_model(ed: EnzymeKatData, template: str = 'inference') -> str:
    paths = {
        k: os.path.join(os.path.dirname(__file__), v)
        for k, v in TEMPLATE_RELATIVE_PATHS.items()
    }
    template_code = read_stan_code_from_path(paths[template])
    functions_block = create_functions_block(ed)
    return functions_block + '\n' + template_code


def create_functions_block(ed: EnzymeKatData) -> str:
    return '\n'.join([
        "functions {",
        "#include big_k_rate_equations.stan",
        "#include haldane_relationships.stan",
        "#include allostery.stan",
        create_fluxes_function(ed),
        create_odes_function(ed),
        create_steady_state_function(),
        "}"
    ])


def create_fluxes_function(ed: EnzymeKatData) -> str:
    mechanism_to_haldane_functions = {
        'uniuni': [create_keq_line],
        'ordered_unibi': [
            create_keq_line,
            create_Kip_ordered_unibi_line,
            create_Kiq_ordered_unibi_line
        ],
    }
    mechanism_to_args_function = {
        'uniuni': get_args_uniuni,
        'ordered_unibi': get_args_ordered_unibi,
        'irr_mass_action': get_args_irr_mass_action,
        'fixed_flux': get_args_fixed_flux
    }
    top = """
    real [] get_fluxes(real[] metabolites, real[] params, real[] known_reals){
    """
    haldane_lines = []
    free_enzyme_ratio_lines = []
    flux_lines = []
    for r in ed.reactions:
        mechanism = r['mechanism']
        if mechanism in mechanism_to_haldane_functions.keys():
            haldane_functions = mechanism_to_haldane_functions[mechanism]
            for f in haldane_functions:
                l = f(ed, r['name'])
                haldane_lines.append(l)
        args_function = mechanism_to_args_function[mechanism]
        args = args_function(ed, r['name'])
        flux_line = f"{mechanism}({args})"
        regulatory_component = create_regulatory_call(ed, r)
        if regulatory_component is not None:
            empty_array_line = "real empty_array[0];"
            if len(free_enzyme_ratio_lines) == 0:
                free_enzyme_ratio_lines.append(empty_array_line)
            free_enzyme_ratio_line = (
                f"real free_enzyme_ratio_{r['name']} = get_free_enzyme_ratio_{mechanism}({args});"
            )
            free_enzyme_ratio_lines.append(free_enzyme_ratio_line)
            flux_line = f"{flux_line} * {regulatory_component}"
        flux_lines.append(flux_line)
    return '\n  '.join([
        "real [] get_fluxes(real[] metabolites, real[] params, real[] known_reals){",
        "\n  ".join(haldane_lines + free_enzyme_ratio_lines),
        "return {",
        ",\n    ".join(flux_lines),
        "};",
        "}",
    ])


def create_odes_function(ed: EnzymeKatData) -> str:
    S = ed.stoichiometry
    fluxes = [f"fluxes[{str(i)}]" for i in range(1, len(S.index) + 1)]
    reaction_to_flux = dict(zip(S.index, fluxes))
    metabolite_lines = {m: '' for m in S.columns}
    for metabolite in S.columns:
        line = metabolite_lines[metabolite]
        for reaction in S.index:
            s = S.loc[reaction, metabolite]
            if s != 0:
                positive_and_not_first = s > 0 and line != ''
                stoich = '+' + str(s) if positive_and_not_first else str(s)
                flux_string = reaction_to_flux[reaction]
                line += f"{stoich}*{flux_string}"
        metabolite_lines[metabolite] = line
    return '\n'.join([
        "real[] get_odes(real[] fluxes){",
        "  return {",
        ",\n    ".join(metabolite_lines.values()),
        "\n  };",
        "}"
    ])
        

def create_steady_state_function():
    return """real[] steady_state_equation(
      real t,
      real[] metabolites,
      real[] params,
      real[] known_reals,
      int[] known_ints
    ){
    for (m in 1:size(metabolites)){
      if (metabolites[m] < 0){
        reject("Metabolite ", m, " is ", metabolites[m], " but should be greater than zero.");
      }
    }
    return get_odes(get_fluxes(metabolites, params, known_reals));
    }
    """
    

def read_stan_code_from_path(path) -> str:
    with open(path, 'r') as f:
        return f.read()


# functions for writing particular lines
def create_keq_line(ed: EnzymeKatData, reaction: str) -> str:
    delta_g_code = ed.parameters.set_index('label').loc[f"{reaction}_delta_g", 'stan_code']
    temperature_code = ed.known_reals.loc['temperature', 'stan_code']
    gas_constant_code = ed.known_reals.loc['gas_constant', 'stan_code']
    return (
        "real {0}_Keq = get_Keq(params[{1}], known_reals[{2}], known_reals[{3}]);"
        .format(reaction, str(delta_g_code), str(temperature_code), str(gas_constant_code))
    )


def create_Kip_ordered_unibi_line(ed: EnzymeKatData, reaction: str) -> str:
    codes = ed.parameters.groupby('label')['stan_code'].first().to_dict()
    kinetic_params = ['Kia', 'Kq', 'Kcat1', 'Kcat2']
    kinetic_param_codes = [codes[reaction + '_' + p] for p in kinetic_params]
    kinetic_params_str = ", ".join([f"params[{str(c)}]" for c in kinetic_param_codes])
    return ''.join([
        f"real {reaction}_Kip = get_Kip_ordered_unibi({reaction}_Keq, ",
        kinetic_params_str,
        ");"
    ])
    

def create_Kiq_ordered_unibi_line(ed: EnzymeKatData, reaction: str) -> str:
    codes = ed.parameters.groupby('label')['stan_code'].first().to_dict()
    kinetic_params = ['Ka', 'Kp', 'Kcat1', 'Kcat2']
    kinetic_param_codes = [codes[reaction + '_' + p] for p in kinetic_params]
    kinetic_params_str = ", ".join([f"params[{str(c)}]" for c in kinetic_param_codes])
    return ''.join([
        f"real {reaction}_Kiq = get_Kiq_ordered_unibi({reaction}_Keq, ",
        kinetic_params_str,
        ");"
    ])


def create_regulatory_call(ed: EnzymeKatData, reaction: dict) -> str:
    param_codes = ed.parameters.groupby('label')['stan_code'].first().to_dict()
    met_codes = ed.metabolite_codes
    if set(reaction.keys()).intersection({'allosteric_inhibitors', 'allosteric_activators'}) == set():
        return None
    else:
        transfer_constant_code = param_codes[reaction['name'] + '_transfer_constant']
        transfer_constant_str = f"params[{transfer_constant_code}]"
        free_enzyme_ratio_str = f"free_enzyme_ratio_{reaction['name']}"
        inhibitors_str = "empty_array"
        dissociation_constant_t_str = "empty_array"
        activators_str = "empty_array"
        dissociation_constant_r_str = "empty_array"
        if 'allosteric_inhibitors' in reaction.keys():
            inhibitor_codes = [met_codes[m] for m in reaction['allosteric_inhibitors']]
            inhibitors_str = "{{{}}}".format(", ".join(map(lambda c: f'metabolites[{str(c)}]', inhibitor_codes)))
            dissociation_constant_t_codes = [
                param_codes[reaction['name'] + '_dissociation_constant_t_' + metabolite]
                for metabolite in reaction['allosteric_inhibitors']
            ]
            dissociation_constant_t_str = "{{{}}}".format(', '.join(map(lambda c: f"params[{c}]", dissociation_constant_t_codes)))
        if 'allosteric_activators' in reaction.keys():
            activator_codes = [met_codes[m] for m in reaction['allosteric_activators']]
            activators_str = "{{{}}}".format(", ".join(map(lambda c: f'metabolites[{str(c)}]', activator_codes)))
            dissociation_constant_t_codes = [
                param_codes[reaction['name'] + '_dissociation_constant_t_' + metabolite]
                for metabolite in reaction['allosteric_activators']
            ]
            dissociation_constant_r_str = "{{{}}}".format(', '.join(map(lambda c: f"params[{c}]", dissociation_constant_r_codes)))
        return ''.join([
            f"get_regulatory_effect(",
            f"{activators_str}, ",
            f"{inhibitors_str}, ",
            f"{free_enzyme_ratio_str}, ",
            f"{dissociation_constant_r_str}, ",
            f"{dissociation_constant_t_str}, ",
            f"{transfer_constant_str}",
            ")"
        ])
    

def get_args_uniuni(ed: EnzymeKatData, reaction: str) -> str:
    kinetic_params = ['Kcat1', 'Kcat2', 'Ka']
    met_codes = ed.metabolite_codes
    param_codes = ed.parameters.groupby('label')['stan_code'].first().to_dict()
    kr_codes = ed.known_reals['stan_code']
    substrate = ed.stoichiometry.loc[reaction].loc[lambda df: df == -1].index[0]
    product = ed.stoichiometry.loc[reaction].loc[lambda df: df == 1].index[0]
    substrate_code, product_code = met_codes[substrate], met_codes[product]
    enzyme_code = kr_codes[f'concentration_{reaction}']
    enzyme_concentration = f"known_reals[{str(enzyme_code)}]"
    kp_codes = {kp: param_codes[f"{reaction}_{kp}"] for kp in kinetic_params}
    kp_strs = {kp: f"params[{str(c)}]" for kp, c in kp_codes.items()}
    for Kcat in ["Kcat1", "Kcat2"]:
        kp_strs[Kcat] = f"{enzyme_concentration}*{kp_strs[Kcat]}"
    S_str = f"metabolites[{str(substrate_code)}]"
    P_str = f"metabolites[{str(product_code)}]"
    Keq_str = f"{reaction}_Keq"
    return ''.join([
        f"{S_str}, {P_str}, {kp_strs['Kcat1']}, ",
        f"{kp_strs['Kcat2']}, {kp_strs['Ka']}, {Keq_str}"
    ]) 


def get_args_ordered_unibi(ed: EnzymeKatData, reaction) -> str:
    kinetic_params = ['Kcat1', 'Kcat2', 'Ka', 'Kp', 'Kq', 'Kia']
    met_codes = ed.metabolite_codes
    param_codes = ed.parameters.groupby('label')['stan_code'].first().to_dict()
    kr_codes = ed.known_reals['stan_code']
    substrate = ed.stoichiometry.loc[reaction].loc[lambda df: df == -1].index[0]
    product_1, product_2 = ed.stoichiometry.loc[reaction].loc[lambda df: df == 1].index
    substrate_code, product_1_code, product_2_code = (
        met_codes[metabolite] for metabolite in [substrate, product_1, product_2]
    )
    enzyme_code = kr_codes[f'concentration_{reaction}']
    enzyme_concentration = f"known_reals[{str(enzyme_code)}]"
    kp_codes = {kp: param_codes[f"{reaction}_{kp}"] for kp in kinetic_params}
    kp_strs = {kp: f"params[{str(c)}]" for kp, c in kp_codes.items()}
    for Kcat in ["Kcat1", "Kcat2"]:
        kp_strs[Kcat] = f"{enzyme_concentration}*{kp_strs[Kcat]}"
    S_str = f"metabolites[{str(substrate_code)}]"
    P1_str = f"metabolites[{str(product_1_code)}]"
    P2_str = f"metabolites[{str(product_2_code)}]"
    Keq_str = f"{reaction}_Keq"
    Kip_str = f"{reaction}_Kip"
    Kiq_str = f"{reaction}_Kiq"
    return ''.join([
        f"{S_str}, {P1_str}, {P2_str}, ",
        f"{kp_strs['Kcat1']}, {kp_strs['Kcat2']}, ",
        f"{kp_strs['Ka']}, {kp_strs['Kp']}, {kp_strs['Kq']}, {kp_strs['Kia']}, ",
        f"{Kip_str}, {Kiq_str}, {Keq_str}",
    ])

def get_args_irr_mass_action(ed: EnzymeKatData, reaction) -> str:
    met_codes = ed.metabolite_codes
    param_codes = ed.parameters.groupby('label')['stan_code'].first().to_dict()
    met = ed.stoichiometry.loc[reaction].loc[lambda df: df != 0].index[0]
    met_code = str(met_codes[met])
    alpha_code = str(param_codes[reaction + '_alpha'])
    met_str = f"metabolites[{met_code}]"
    alpha_str = f"params[{alpha_code}]"
    return f"{met_str}, {alpha_str}"

def get_args_fixed_flux(ed: EnzymeKatData, reaction) -> str:
    kr_code = ed.known_reals.loc[f'concentration_{reaction}', 'stan_code']
    kr_str = f"known_reals[{str(kr_code)}]"
    return f"{kr_str}"
