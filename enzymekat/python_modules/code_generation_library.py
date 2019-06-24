import pandas as pd
from python_modules.enzymekat_data import EnzymeKatData

# functions for writing particular lines
def create_keq_line(ed: EnzymeKatData, reaction: str) -> str:
    delta_g_code = (
        ed.parameters
        .set_index(['reaction', 'parameter'])
        .loc[(reaction, 'delta_g'), 'stan_code']
        .astype(str)
    )
    temperature_code, gas_constant_code = (
        ed.known_reals.loc[['temperature', 'gas_constant'], 'stan_code']
        .astype(str)
    )
    return (
        "real {0}_Keq = get_Keq(params[{1}], params[{2}], params[{3}]);"
        .format(reaction, delta_g_code, temperature_code, gas_constant_code)
    )


def create_Kip_ordered_unibi_line(ed: EnzymeKatData, reaction: str) -> str:
    params = ['Kia', 'Kq', 'Kcat1', 'Kcat2']
    param_codes = (
        ed.parameters
        .set_index(['reaction', 'parameter'])
        .loc[reaction].loc[params, 'stan_code']
        .astype(str)
        .tolist()
    )
    params_component = ", ".join([f"params[{c}]" for c in param_codes])
    return (
        f"real {reaction}_Kip = get_Kip_ordered_unibi({reaction}_Keq, "
        + params_component
        + ");"
    )
    

def create_Kiq_ordered_unibi_line(ed: EnzymeKatData, reaction: str) -> str:
    params = ['Ka', 'Kp', 'Kcat1', 'Kcat2']
    param_codes = (
        ed.parameters
        .set_index(['reaction', 'parameter'])
        .loc[reaction].loc[params, 'stan_code']
        .astype(str)
        .tolist()
    )
    params_component = ", ".join([f"params[{c}]" for c in param_codes])
    return (
        f"real {reaction}_Kiq = get_Kiq_ordered_unibi({reaction}_Keq, "
        + params_component
        + ");"
    )


def create_uniuni_call(ed: EnzymeKatData, reaction) -> str:
    substrate = ed.stoichiometry.loc[reaction].loc[lambda df: df == -1].index[0]
    product = ed.stoichiometry.loc[reaction].loc[lambda df: df == 1].index[0]
    substrate_code, product_code = (
        ed.metabolites
        .set_index('name')
        .loc[[substrate, product], 'stan_code']
        .astype(str)
    )
    enzyme_code = ed.known_reals.loc[f'concentration_{reaction}', 'stan_code']
    Kcat1_code, Kcat2_code, Ka_code = (
        ed.parameters.set_index(['reaction', 'parameter'])
        .loc[reaction].loc[['Kcat1', 'Kcat2', 'Ka'], 'stan_code']
        .astype(str)
    )
    return (
        f"uniuni(metabolites[{substrate_code}], metabolites[{product_code}], "
        + f"known_reals[{enzyme_code}]*params[{Kcat1_code}], "
        + f"known_reals[{enzyme_code}]*params[{Kcat2_code}], "
        + f"params[{Ka_code}], "
        + f"{reaction}_Keq);"
    )


def create_ordered_unibi_call(ed: EnzymeKatData, reaction) -> str:
    non_multiplied_params = ['Ka', 'Kp', 'Kq', 'Kia']
    haldane_params = ['Kip', 'Kiq', 'Keq']
    substrate = ed.stoichiometry.loc[reaction].loc[lambda df: df == -1].index[0]
    product_1, product_2 = ed.stoichiometry.loc[reaction].loc[lambda df: df == 1].index
    substrate_code, product_1_code, product_2_code = (
        ed.metabolites
        .set_index('name')
        .loc[[substrate, product_1, product_2], 'stan_code']
        .astype(str)
    )
    enzyme_code = ed.known_reals.loc[f'concentration_{reaction}', 'stan_code']
    Kcat1_code, Kcat2_code = (
        ed.parameters.set_index(['reaction', 'parameter'])
        .loc[reaction].loc[['Kcat1', 'Kcat2'], 'stan_code']
        .astype(str)
    )
    non_multiplied_param_codes  = (
        ed.parameters.set_index(['reaction', 'parameter'])
        .loc[reaction].loc[non_multiplied_params, 'stan_code']
        .astype(str)
    )
    return (
        "ordered_unibi("
        + f"metabolites[{substrate_code}], metabolites[{product_1_code}], metabolites[{product_2_code}], "
        + f"known_reals[{enzyme_code}]*params[{Kcat1_code}], "
        + f"known_reals[{enzyme_code}]*params[{Kcat2_code}], "
        + ", ".join([f"params[c]" for c in non_multiplied_param_codes])
        + ", "
        + ", ".join([f"{reaction}_{p}" for p in haldane_params])
        + ");"
    )

def create_irr_mass_action_call(ed: EnzymeKatData, reaction) -> str:
    metabolite = ed.stoichiometry.loc[reaction].loc[lambda df: df != 1].index[0]
    metabolite_code = (
        ed.metabolites
        .set_index('name')
        .loc[metabolite, 'stan_code']
        .astype(int)
    )
    alpha_code = (
        ed.parameters
        .set_index(['reaction', 'parameter'])
        .loc[(reaction, 'alpha'), 'stan_code']
        .astype(str)
    )
    return (
        f"irr_mass_action(metabolites[{metabolite_code}], params[{alpha_code}]);"
    )





