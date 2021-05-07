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

"""General purpose utility functions."""

from typing import Dict, Iterable
import collections

import numpy as np
import pandas as pd
import sympy as sp
from scipy.stats import norm


def codify(lx: Iterable[str]) -> Dict[str, int]:
    """Turn an iterable of strings into a dictionary mapping them to integer indexes."""
    return dict(zip(lx, range(1, len(lx) + 1)))


def get_lognormal_parameters_from_quantiles(x1, p1, x2, p2):
    """Find parameters for a lognormal distribution from two quantiles.

    i.e. get mu and sigma such that if X ~ lognormal(mu, sigma), then pr(X <
    x1) = p1 and pr(X < x2) = p2.

    """
    logx1 = np.log(x1)
    logx2 = np.log(x2)
    denom = norm.ppf(p2) - norm.ppf(p1)
    sigma = (logx2 - logx1) / denom
    mu = (logx1 * norm.ppf(p2) - logx2 * norm.ppf(p1)) / denom
    return mu, sigma


def get_normal_parameters_from_quantiles(x1, p1, x2, p2):
    """Find parameters for a normal distribution from two quantiles.

    i.e. get mu and sigma such that if X ~ normal(mu, sigma), then pr(X <
    x1) = p1 and pr(X < x2) = p2.

    """
    denom = norm.ppf(p2) - norm.ppf(p1)
    sigma = (x2 - x1) / denom
    mu = (x1 * norm.ppf(p2) - x2 * norm.ppf(p1)) / denom
    return mu, sigma


def get_null_space(a, rtol=1e-5):
    """Calulate the null space of a matrix."""
    u, s, v = np.linalg.svd(a)
    rank = (s > rtol * s[0]).sum()
    return v[rank:].T.copy()


def get_rref(mat):
    """Return reduced row echelon form of a matrix."""
    return sp.Matrix(mat).rref(iszerofunc=lambda x: abs(x) < 1e-10)[0]


def export_scaled_params_from_draws(infd, chain, draw):
    """Extact log centered parameters from an infd object."""
    list_of_input_inits = [
        "log_km_z",
        "drain_z",
        "log_ki_z",
        "log_dissociation_constant_t_z",
        "log_dissociation_constant_r_z",
        "log_transfer_constant_z",
        "log_kcat_z",
        "log_phos_kcat_z",
        "log_conc_unbalanced_z",
        "log_enzyme_z",
        "log_phos_conc_z",
        "formation_energy_z",
    ]

    input_dict = {
        par_name: infd.posterior[par_name][chain, draw].values.tolist()
        for par_name in list_of_input_inits
        if par_name in infd.posterior.variables.keys()
    }
    return input_dict


def export_params_from_draw(infd, chain, draw):
    """Extact parameters from an infd object."""
    list_of_input_inits = [
        "km",
        "drain",
        "ki",
        "diss_t",
        "diss_r",
        "transfer_constant",
        "kcat",
        "phos_enzyme_kcat",
        "conc_unbalanced",
        "enzyme",
        "phos_enzyme_conc",
        "formation_energy",
    ]

    input_dict = {
        par_name: infd.posterior[par_name][chain, draw].values.tolist()
        for par_name in list_of_input_inits
        if par_name in infd.posterior.variables.keys()
    }
    return input_dict

def return_ordered_combinations(list1, list2):
    ordered_combinations = []
    for el1 in list1:
        for el2 in list2:
            ordered_combinations.append([el1, el2])
    if len(list2) == 0:
        ordered_combinations = [[], []]
    return np.transpose(ordered_combinations)

def collect_input_coords(scs):
    output_df = pd.DataFrame()

    parameter_names_and_coords = {
        "kcat": [["enzyme_id"], [scs.enzymes]],
        "km": [["enzyme_id", "mic_id"], [scs.km_enzs, scs.km_mics]],
        "formation_energy": [["metabolite_id"], [scs.metabolites]],
        "diss_t": [["enzyme_id", "mic_id"], [scs.ai_enzs, scs.ai_mics]],
        "diss_r": [["enzyme_id", "mic_id"], [scs.aa_enzs, scs.aa_mics]],
        "ki": [["enzyme_id", "mic_id"], [scs.ci_enzs, scs.ci_mics]],
        "transfer_constant": [["enzyme_id"], [scs.allosteric_enzymes]],
        "phos_enzyme_kcat": [["phos_enz_id"], [scs.phos_enzs]],
        "conc_unbalanced": [["experiment_id", "mic_id"], return_ordered_combinations(scs.experiments, scs.unbalanced_mics)],
        "enzyme": [["experiment_id", "enzyme_id"], return_ordered_combinations(scs.experiments, scs.enzymes)],
        "drain": [["experiment_id", "drain_id"], return_ordered_combinations(scs.experiments, scs.drains)],
        "phos_enzyme_conc": [["experiment_id", "phos_enz_id"], return_ordered_combinations(scs.experiments, scs.phos_enzs)],
    }

    for par, idx in parameter_names_and_coords.items():
        if len(idx[1]) == 1:
            tmp_df = pd.DataFrame.from_dict({"parameter_name": par,
                                             idx[0][0]: idx[1][0]})
            output_df = output_df.append(tmp_df, ignore_index=True)
        if len(idx[1]) == 2:
            tmp_df = pd.DataFrame.from_dict({"parameter_name": par,
                                             idx[0][0]: idx[1][0],
                                             idx[0][1]: idx[1][1]})
            output_df = output_df.append(tmp_df)

    return output_df


def convert_normal_to_scaled_values(input_dict, priors, stan_coords):
    """Converts values on normal scale to the scaled values.
    
    :input_dict: dictionary of parameter values.
    :priors: dictionary of priors.
    """
    Init_variable = collections.namedtuple("Init_variable", "in_id out_id p norm_scale exp_dep")
    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', 20)
    input_coords = collect_input_coords(stan_coords)
    coord_names = ["parameter_name", "experiment_id", "mic_id", "enzyme_id", "drain_id", "metabolite_id"]

    init_df = input_coords.merge(input_dict,
                                   left_on=coord_names,
                                   right_on=coord_names,
                                   how="inner")
    
    list_of_input_inits = [
        Init_variable(in_id="km", out_id = "log_km_z", p=priors.km_priors.location, norm_scale=False, exp_dep=False),
        Init_variable(in_id="drain", out_id = "drain_z", p=priors.drain_priors.location, norm_scale=True, exp_dep=True),
        Init_variable(in_id="ki", out_id = "log_ki_z", p=priors.inhibition_constant_priors.location, norm_scale=False, exp_dep=False),
        Init_variable(in_id="diss_t", out_id = "log_dt_z", p=priors.tense_dissociation_constant_priors.location, norm_scale=False, exp_dep=False),
        Init_variable(in_id="diss_r", out_id = "log_dr_z", p=priors.relaxed_dissociation_constant_priors.location, norm_scale=False, exp_dep=False),
        Init_variable(in_id="transfer_constant", out_id = "log_tc_z", p=priors.transfer_constant_priors.location, norm_scale=False, exp_dep=False),
        Init_variable(in_id="kcat", out_id = "log_kcat_z", p=priors.kcat_priors.location, norm_scale=False, exp_dep=False),
        Init_variable(in_id="phos_enzyme_kcat", out_id = "log_phos_kcat_z", p=priors.phos_kcat_priors.location, norm_scale=False, exp_dep=False),
        Init_variable(in_id="conc_unbalanced", out_id = "log_conc_unbalanced_z", p=priors.unbalanced_metabolite_priors.location, norm_scale=False, exp_dep=True),
        Init_variable(in_id="enzyme", out_id = "log_enzyme_z", p=priors.enzyme_concentration_priors.location, norm_scale=False, exp_dep=True),
        Init_variable(in_id="phos_enzyme_conc", out_id = "log_phos_conc_z", p=priors.phos_enz_concentration_priors.location, norm_scale=False, exp_dep=True),
        Init_variable(in_id="formation_energy", out_id = "fe_z", p=priors.formation_energy_priors.location, norm_scale=True, exp_dep=False),
    ]

    inits = {}

    for init in list_of_input_inits:
        if init.in_id in input_dict["parameter_name"].unique():
            if init.exp_dep:
                prior_location = init.p.stack().reset_index().rename(columns={0: "value"})
                prior_location["location"] = prior_location["value"]
                tmp_df = input_coords.loc[input_coords["parameter_name"] == init.in_id]
                parameter_values = prior_location.merge(tmp_df,
                                                  left_on=prior_location.columns.values.tolist()[:-2],
                                                  right_on=prior_location.columns.values.tolist()[:-2])
                print(parameter_values)
                init_input = parameter_values.merge(init_df,
                                                    left_on=coord_names,
                                                    right_on=coord_names)
                print(init_input)

                input_array = [[parameter_values.loc[parameter_values["experiment_id"] == exp]["value"]] for exp in stan_coords.experiments]
                prior_array = [[parameter_values.loc[parameter_values["experiment_id"] == exp]["location"]] for exp in stan_coords.experiments]
                if init.norm_scale:
                    inits[init.out_id] = z_for_mat(input_array, prior_array)
                else:
                    inits[init.out_id] = logz_for_mat(input_array, prior_array)
            # else:
            #     prior_location = init.p.reset_index()
            #     tmp_df = input_coords.loc[input_coords["parameter_name"] == init.in_id]
            #     parameter_values = prior_location.merge(tmp_df,
            #                                       left_on=prior_location.columns.values.tolist()[:-1],
            #                                       right_on=prior_location.columns.values.tolist()[:-1])
            #     if init.norm_scale:
            #         inits[init.out_id] = z_for_vec(parameter_values["value"], parameter_values["location"])
            #     else:
            #         inits[init.out_id] = logz_for_vec(parameter_values["value"], parameter_values["location"])

    return inits

def logz_for_vec(input_value, priors):
    if len(input_value) == 0:
        return []
    iv = np.array(input_value)
    lc = np.array(priors[0])
    s = np.array(priors[1])
    return (np.log(iv) - np.log(lc)) / s

def z_for_vec(input_value, priors):
    if len(input_value) == 0:
        return []
    iv = np.array(input_value)
    lc = np.array(priors[0])
    s = np.array(priors[1])
    return (iv - lc) / s

def logz_for_mat(input_value, priors):
    if len(input_value) == 0:
        return np.array([[] for ex in priors])
    out = []
    for i, [lc_v, s_v] in enumerate(zip(priors[0], priors[1])):
        iv = np.array(input_value[i])
        lc = np.array(lc_v)
        s = np.array(s_v)
        out.append((np.log(iv) - np.log(lc)) / s)
    return np.array(out)

def z_for_mat(input_value, priors):
    if len(input_value) == 0:
        return np.array([[] for ex in priors])
    out = []
    for i, [lc_v, s_v] in enumerate(zip(priors[0], priors[1])):
        iv = np.array(input_value[i])
        lc = np.array(lc_v)
        s = np.array(s_v)
        out.append((iv - lc) / s)
    return np.array(out)
