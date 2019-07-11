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
    return (
        "functions {\n"
        + "#include big_k_rate_equations.stan\n"
        + "#include haldane_relationships.stan\n"
        + create_fluxes_function(ed)
        + "\n"
        + create_odes_function(ed)
        + "\n"
        + create_steady_state_function()
        + "}\n"
    )


def create_fluxes_function(ed: EnzymeKatData) -> str:
    mechanism_to_haldane_functions = {
        'uniuni': [create_keq_line],
        'ordered_unibi': [
            create_keq_line,
            create_Kip_ordered_unibi_line,
            create_Kiq_ordered_unibi_line
        ],
    }
    mechanism_to_flux_function = {
        'uniuni': create_uniuni_call,
        'ordered_unibi': create_ordered_unibi_call,
        'irr_mass_action': create_irr_mass_action_call,
        'fixed_flux': create_fixed_flux_call
    }

    top = """
    real [] get_fluxes(real[] metabolites, real[] params, real[] known_reals){
    """
    haldane_lines = []
    flux_lines = []
    for r in ed.reactions:
        mechanism = r['mechanism']
        if mechanism in mechanism_to_haldane_functions.keys():
            haldane_functions = mechanism_to_haldane_functions[mechanism]
            for f in haldane_functions:
                l = f(ed, r['name'])
                haldane_lines.append(l)
        flux_line = mechanism_to_flux_function[mechanism](ed, r['name'])
        flux_lines.append(flux_line)
    return (
        "real [] get_fluxes(real[] metabolites, real[] params, real[] known_reals){\n  "
        + "\n  ".join(haldane_lines)
        + "\n  return {\n    "
        + ",\n    ".join(flux_lines)
        + "\n  };"
        + "\n}"
    )
    

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
    return (
        "real[] get_odes(real[] fluxes){\n"
        + "  return {\n    "
        + ",\n    ".join(metabolite_lines.values())
        + "\n  };\n"
        + "}"
    )
        

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
    delta_g_code = str(ed.stan_codes['thermodynamic_parameter'][reaction + '_delta_g'])
    temperature_code = str(ed.stan_codes['known_real']['temperature'])
    gas_constant_code = str(ed.stan_codes['known_real']['gas_constant'])
    return (
        "real {0}_Keq = get_Keq(params[{1}], known_reals[{2}], known_reals[{3}]);"
        .format(reaction, delta_g_code, temperature_code, gas_constant_code)
    )


def create_Kip_ordered_unibi_line(ed: EnzymeKatData, reaction: str) -> str:
    kinetic_params = ['Kia', 'Kq', 'Kcat1', 'Kcat2']
    N_therm = len(ed.thermodynamic_parameters)
    kinetic_param_codes = [
        str(N_therm + ed.stan_codes['kinetic_parameter'][reaction + '_' + param])
        for param in kinetic_params
    ]
    kinetic_params_component = ", ".join([f"params[{c}]" for c in kinetic_param_codes])
    return (
        f"real {reaction}_Kip = get_Kip_ordered_unibi({reaction}_Keq, "
        + kinetic_params_component
        + ");"
    )
    

def create_Kiq_ordered_unibi_line(ed: EnzymeKatData, reaction: str) -> str:
    params = ['Ka', 'Kp', 'Kcat1', 'Kcat2']
    N_therm = len(ed.thermodynamic_parameters)
    kinetic_param_codes = [
        str(N_therm + ed.stan_codes['kinetic_parameter'][reaction + '_' + param])
        for param in params
    ]
    kinetic_params_component = ", ".join([f"params[{c}]" for c in kinetic_param_codes])
    return (
        f"real {reaction}_Kiq = get_Kiq_ordered_unibi({reaction}_Keq, "
        + kinetic_params_component
        + ");"
    )


def create_uniuni_call(ed: EnzymeKatData, reaction) -> str:
    N_therm = len(ed.thermodynamic_parameters)
    substrate = ed.stoichiometry.loc[reaction].loc[lambda df: df == -1].index[0]
    product = ed.stoichiometry.loc[reaction].loc[lambda df: df == 1].index[0]
    kinetic_params = ['Kcat1', 'Kcat2', 'Ka']
    substrate_code, product_code = (
        str(ed.stan_codes['metabolite'][metabolite])
        for metabolite in [substrate, product]
    )
    enzyme_code = ed.stan_codes['known_real'][f'concentration_{reaction}']
    Kcat1_code, Kcat2_code, Ka_code = [
        str(N_therm + ed.stan_codes['kinetic_parameter'][reaction + '_' + param])
        for param in kinetic_params
    ]
    return (
        f"uniuni(metabolites[{substrate_code}], metabolites[{product_code}], "
        + f"known_reals[{enzyme_code}]*params[{Kcat1_code}], "
        + f"known_reals[{enzyme_code}]*params[{Kcat2_code}], "
        + f"params[{Ka_code}], "
        + f"{reaction}_Keq)"
    )


def create_ordered_unibi_call(ed: EnzymeKatData, reaction) -> str:
    N_therm = len(ed.thermodynamic_parameters)
    non_multiplied_kinetic_params = ['Ka', 'Kp', 'Kq', 'Kia']
    haldane_params = ['Kip', 'Kiq', 'Keq']
    substrate = ed.stoichiometry.loc[reaction].loc[lambda df: df == -1].index[0]
    product_1, product_2 = ed.stoichiometry.loc[reaction].loc[lambda df: df == 1].index
    substrate_code, product_1_code, product_2_code = [
        str(ed.stan_codes['metabolite'][metabolite])
        for metabolite in [substrate, product_1, product_2]
    ]
    enzyme_code = ed.stan_codes['known_real'][f'concentration_{reaction}']
    Kcat1_code, Kcat2_code = (
        str(N_therm + ed.stan_codes['kinetic_parameter'][reaction + '_' + param])
        for param in ['Kcat1', 'Kcat2']
    )
    non_multiplied_kinetic_param_codes  = (
        str(N_therm + ed.stan_codes['kinetic_parameter'][reaction + '_' + param])
        for param in non_multiplied_kinetic_params
    )
    return (
        "ordered_unibi("
        + f"metabolites[{substrate_code}], metabolites[{product_1_code}], metabolites[{product_2_code}], "
        + f"known_reals[{enzyme_code}]*params[{Kcat1_code}], "
        + f"known_reals[{enzyme_code}]*params[{Kcat2_code}], "
        + ", ".join([f"params[{c}]" for c in non_multiplied_kinetic_param_codes])
        + ", "
        + ", ".join([f"{reaction}_{p}" for p in haldane_params])
        + ")"
    )

def create_irr_mass_action_call(ed: EnzymeKatData, reaction) -> str:
    N_therm = len(ed.thermodynamic_parameters)
    metabolite = ed.stoichiometry.loc[reaction].loc[lambda df: df != 0].index[0]
    metabolite_code = str(ed.stan_codes['metabolite'][metabolite])
    alpha_code = str(N_therm + ed.stan_codes['kinetic_parameter'][reaction + '_alpha'])
    return (
        f"irr_mass_action(metabolites[{metabolite_code}], params[{alpha_code}])"
    )

def create_fixed_flux_call(ed: EnzymeKatData, reaction) -> str:
    N_therm = len(ed.thermodynamic_parameters)
    param_code = str(N_therm + ed.stan_codes['kinetic_parameter'][reaction + '_v'])
    return (
        f"fixed_flux(params[{param_code}])"
    )
