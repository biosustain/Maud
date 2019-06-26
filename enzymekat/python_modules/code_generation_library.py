import pandas as pd
from python_modules.enzymekat_data import EnzymeKatData

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
