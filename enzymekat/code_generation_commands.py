import os
import pandas as pd
import code_generation_library as lib
from enzymekat_data import EnzymeKatData


RELATIVE_PATHS = {
    'inference_model_template': 'stan_code/inference_model_template.stan',
    'steady_state_equation': 'stan_code/steady_state_equation.stan'
}
PATHS = {k: os.path.join(os.path.dirname(__file__), v) for k, v in RELATIVE_PATHS.items()}
MECHANISM_TO_HALDANE_FUNCTIONS = {
    'uniuni': [lib.create_keq_line],
    'ordered_unibi': [
        lib.create_keq_line,
        lib.create_Kip_ordered_unibi_line,
        lib.create_Kiq_ordered_unibi_line
    ],
}
MECHANISM_TO_FLUX_FUNCTION = {
    'uniuni': lib.create_uniuni_call,
    'ordered_unibi': lib.create_ordered_unibi_call,
    'irr_mass_action': lib.create_irr_mass_action_call,
    'fixed_flux': lib.create_fixed_flux_call
}


def create_stan_model(ed: EnzymeKatData) -> str:
    functions_block = create_functions_block(ed)
    inference_model_template = read_stan_code_from_path(PATHS['inference_model_template'])
    return functions_block + '\n' + inference_model_template


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
    top = """
    real [] get_fluxes(real[] metabolites, real[] params, real[] known_reals){
    """
    haldane_lines = []
    flux_lines = []
    for r in ed.reactions:
        mechanism = r['mechanism']
        if mechanism in MECHANISM_TO_HALDANE_FUNCTIONS.keys():
            haldane_functions = MECHANISM_TO_HALDANE_FUNCTIONS[mechanism]
            for f in haldane_functions:
                l = f(ed, r['name'])
                haldane_lines.append(l)
        flux_line = MECHANISM_TO_FLUX_FUNCTION[mechanism](ed, r['name'])
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
