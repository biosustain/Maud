import os
import warnings
from typing import List, Union

import cmdstanpy
import numpy as np
import pandas as pd

from maud.data_model import (
    IndPrior1d,
    IndPrior2d,
    MaudInput,
    MeasurementSet,
    PriorSet,
    StanCoordSet,
)
from maud.utils import codify, get_null_space, get_rref, extract_inits_from_infd
from maud.analysis import load_infd
from maud.sampling import get_stoics, get_phos_act_inh_matrix, \
                          _get_km_lookup, _get_conc_init, get_config_dict, \
                          validate_specified_fluxes

HERE = os.path.dirname(os.path.abspath(__file__))
INPUT_NAME = "../../methionine_input"
INPUT_FOLDER = os.path.join(HERE, INPUT_NAME)
INPUT_CSV = os.path.join(INPUT_FOLDER, "simulation_input_values.csv")
INCLUDE_PATH = ""
DEFAULT_PRIOR_LOC_DRAIN = None
DEFAULT_PRIOR_SCALE_DRAIN = None
STAN_PROGRAM_RELATIVE_PATH = "simulate_ensemble.stan"

DEFAULT_ODE_CONFIG = {
    "rel_tol": 1e-9,
    "abs_tol": 1e-9,
    "max_num_steps": int(1e9),
    "timepoint": 500,
}

SIM_CONFIG = {
    "chains": 1,
    "fixed_param": True,
    "inits": 0,
    "iter_warmup": 0,
    "iter_sampling": 1,
    "show_progress": False,
    "threads_per_chain": 1,
}


def simulate(mi: MaudInput, output_dir: str):
    """Generate simulations from the prior mean.

    :param mi: a MaudInput object
    :param output_dir: a string specifying where to save the output.
    """
    config = {**SIM_CONFIG, **{"output_dir": output_dir}}
    return _simulate_given_input(mi, output_dir, config)

def get_simulation_input_data(mi: MaudInput, csv) -> dict:
    """Get the input to inference_model.stan from a MaudInput object.

    :param mi: a MaudInput object

    """
    validate_specified_fluxes(mi)
    sorted_enzymes = sorted(
        [e for r in mi.kinetic_model.reactions for e in r.enzymes],
        key=lambda e: codify(mi.stan_coords.enzymes)[e.id],
    )
    sorted_mics = sorted(
        mi.kinetic_model.mics,
        key=lambda m: codify(mi.stan_coords.mics)[m.id],
    )

    water_stoichiometry = [e.water_stoichiometry for e in sorted_enzymes]
    mic_to_met = [
        codify(mi.stan_coords.metabolites)[mic.metabolite_id] for mic in sorted_mics
    ]
    S_enz, S_to_flux, S_full, S_drain, _ = get_stoics(mi)
    knockout_matrix_enzyme = mi.measurements.enz_knockouts.astype(int)
    knockout_matrix_phos = mi.measurements.phos_knockouts.astype(int)
    S_phos_act, S_phos_inh = get_phos_act_inh_matrix(mi)
    unbalanced_mic_ix, balanced_mic_ix = (
        [codify(mi.stan_coords.mics)[mic_id] for mic_id in ix]
        for ix in (
            mi.stan_coords.unbalanced_mics,
            mi.stan_coords.balanced_mics,
        )
    )
    allosteric_enzymes = [e for e in sorted_enzymes if e.allosteric]
    config_dict = get_config_dict(mi)
    mic_ix = codify(mi.stan_coords.mics)
    return {
        **{
            # sizes
            "N_samples": 1,
            "N_mic": len(mi.kinetic_model.mics),
            "N_unbalanced": len(mi.stan_coords.unbalanced_mics),
            "N_metabolite": len(mi.stan_coords.metabolites),
            "N_km": len(mi.stan_coords.km_mics),
            "N_reaction": len(mi.stan_coords.reactions),
            "N_enzyme": len(mi.stan_coords.enzymes),
            "N_phosphorylation_enzymes": len(mi.stan_coords.phos_enzs),
            "N_experiment": len(mi.stan_coords.experiments),
            "N_flux_measurement": len(mi.stan_coords.yflux_rxns),
            "N_enzyme_measurement": len(mi.stan_coords.yenz_enzs),
            "N_conc_measurement": len(mi.stan_coords.yconc_mics),
            "N_ki": len(mi.stan_coords.ci_mics),
            "N_ai": len(mi.stan_coords.ai_mics),
            "N_aa": len(mi.stan_coords.aa_mics),
            "N_ae": len(allosteric_enzymes),
            "N_drain": len(mi.stan_coords.drains),
            # indexes
            "unbalanced_mic_ix": unbalanced_mic_ix,
            "balanced_mic_ix": balanced_mic_ix,
            "ci_ix": [mic_ix[m] for m in mi.stan_coords.ci_mics],
            "ai_ix": [mic_ix[m] for m in mi.stan_coords.ai_mics],
            "aa_ix": [mic_ix[m] for m in mi.stan_coords.aa_mics],
            # network properties
            "S_enz": S_enz.T.values,
            "S_to_flux_map": S_to_flux.values,
            "S_drain": S_drain.T.values,
            "S_full": S_full.T.values,
            "S_phos_act": S_phos_act,
            "S_phos_inh": S_phos_inh,
            "water_stoichiometry": water_stoichiometry,
            "mic_to_met": mic_to_met,
            "km_lookup": _get_km_lookup(mi),
            "is_knockout": knockout_matrix_enzyme.values.tolist(),
            "is_phos_knockout": knockout_matrix_phos.values.tolist(),
            "subunits": [e.subunits for e in sorted_enzymes],
            "n_ci": [len(e.modifiers["competitive_inhibitor"]) for e in sorted_enzymes],
            "n_ai": [len(e.modifiers["allosteric_inhibitor"]) for e in sorted_enzymes],
            "n_aa": [len(e.modifiers["allosteric_activator"]) for e in sorted_enzymes],
        },
        **extract_inits_from_infd(mi, csv),
        **config_dict,
    }

def _simulate_given_input(mi: MaudInput, output_dir: str, config: dict):
    """Call CmdStanModel.sample, having already specified all arguments.

    :param mi: a MaudInput object
    :param output_dir: a string specifying where to save the output.
    :param config: a dictionary of keyword arguments to CmdStanModel.sample.
    """

    input_filepath = os.path.join(output_dir, "input_data.json")
    input_data = get_simulation_input_data(mi, INPUT_CSV)
    cmdstanpy.utils.jsondump(input_filepath, input_data)
    stan_program_filepath = os.path.join(HERE, STAN_PROGRAM_RELATIVE_PATH)
    include_path = os.path.join(HERE, INCLUDE_PATH)
    cpp_options = {}
    stanc_options = {"include_paths": [include_path]}
    if config["threads_per_chain"] != 1:
        cpp_options["STAN_THREADS"] = True
        os.environ["STAN_NUM_THREADS"] = str(config["threads_per_chain"])
    model = cmdstanpy.CmdStanModel(
        stan_file=stan_program_filepath,
        stanc_options=stanc_options,
        cpp_options=cpp_options,
    )
    stanfit = model.sample(data=input_filepath, **config)
    infd = load_infd(stanfit, mi)
    infd.to_netcdf(os.path.join(output_dir, "..", "infd.ncdf"))
    return stanfit, infd