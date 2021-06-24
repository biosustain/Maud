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

"""Code for sampling from a posterior distribution."""

import itertools
import os
import warnings
from typing import Dict, List, Union

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
from maud.utils import codify, get_null_space, get_rref


HERE = os.path.dirname(os.path.abspath(__file__))
INCLUDE_PATH = ""
DEFAULT_PRIOR_LOC_DRAIN = None
DEFAULT_PRIOR_SCALE_DRAIN = None
STAN_PROGRAM_RELATIVE_PATH = "model.stan"

DEFAULT_SAMPLE_CONFIG = {
    "iter_warmup": 5,
    "iter_sampling": 5,
    "chains": 2,
    "max_treedepth": 11,
    "show_progress": True,
    "step_size": 0.025,
    "adapt_delta": 0.99,
    "save_warmup": True,
    "threads_per_chain": 1,
}

SIM_CONFIG = {
    "chains": 1,
    "fixed_param": True,
    "iter_warmup": 0,
    "show_progress": False,
    "threads_per_chain": 1,
}


def sample(mi: MaudInput, output_dir: str) -> cmdstanpy.CmdStanMCMC:
    """Sample from the posterior defined by mi.

    :param mi: a MaudInput object
    :param output_dir: a string specifying where to save the output.
    """
    config = {
        **DEFAULT_SAMPLE_CONFIG,
        **mi.config.cmdstanpy_config,
        **{"output_dir": output_dir},
    }
    return _sample_given_config(mi, output_dir, config)


def simulate(mi: MaudInput, output_dir: str, n: int) -> cmdstanpy.CmdStanMCMC:
    """Generate simulations from the prior mean.

    :param mi: a MaudInput object
    :param output_dir: a string specifying where to save the output.
    """
    config = {**SIM_CONFIG, **{"output_dir": output_dir, "iter_sampling": n}}
    return _sample_given_config(mi, output_dir, config)


def _sample_given_config(
    mi: MaudInput, output_dir: str, config: dict
) -> cmdstanpy.CmdStanMCMC:
    """Call CmdStanModel.sample, having already specified all arguments.

    :param mi: a MaudInput object
    :param output_dir: a string specifying where to save the output.
    :param config: a dictionary of keyword arguments to CmdStanModel.sample.
    """

    input_filepath = os.path.join(output_dir, "input_data.json")
    inits_filepath = os.path.join(output_dir, "inits.json")
    coords_filepath = os.path.join(output_dir, "coords.json")
    input_data = get_input_data(mi)
    inits = {k: v.values for k, v in mi.inits.items()}
    cmdstanpy.utils.jsondump(input_filepath, input_data)
    cmdstanpy.utils.jsondump(inits_filepath, inits)
    cmdstanpy.utils.jsondump(coords_filepath, mi.stan_coords.__dict__)
    config["inits"] = inits_filepath
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
    return model.sample(data=input_filepath, **config)


def get_stoichiometry(mi: MaudInput) -> pd.DataFrame:
    """Get a stoichiometric matrix dataframe from a MaudInput.

    :param mi: A MaudInput object
    """
    S = pd.DataFrame(0, index=mi.stan_coords.edges, columns=mi.stan_coords.mics)
    for rxn in mi.kinetic_model.reactions:
        if rxn.reaction_mechanism == "drain":
            for met, stoic in rxn.stoichiometry.items():
                S.loc[rxn.id, met] = stoic
        else:
            for enz in rxn.enzymes:
                for met, stoic in rxn.stoichiometry.items():
                    S.loc[enz.id, met] = stoic
    return S.T


def validate_specified_fluxes(mi: MaudInput):
    """Check that appropriate fluxes have been measured.

    :param mi: A MaudInput object
    """
    S = get_stoichiometry(mi)
    balanced_mic_ix = mi.stan_coords.balanced_mics
    reactions = mi.stan_coords.reactions
    for exp in mi.stan_coords.experiments:
        measured_rxns = mi.measurements.yflux.reset_index()["target_id"].unique()
        flux_paths = get_null_space(S.loc[balanced_mic_ix].values)
        _, n_dof = np.shape(flux_paths)
        rref_flux_paths = np.matrix(get_rref(flux_paths.T))
        rref_flux_paths[np.abs(rref_flux_paths) < 1e-10] = 0
        flux_paths_df = pd.DataFrame(rref_flux_paths, columns=reactions)
        for _, flux_path in flux_paths_df.iterrows():
            if any(flux_paths[measured_rxns]) != 0:
                pass
            else:
                possible_measurements = []
                for rxn, st in flux_path.items():
                    if st != 0:
                        possible_measurements.append(rxn)
                msg = (
                    "\nYour system appears to be underdetermined in "
                    + f"experiment: {exp}\n"
                    + "Please define a reaction from the following list:\n"
                    + "\n".join(possible_measurements)
                )
                warnings.warn(msg)

        if len(measured_rxns) > n_dof:
            msg = (
                "You appear to have specified too many reactions.\n"
                + "This will bias the statistical model\n"
                + "as the measurements are not independent."
            )
            warnings.warn(msg)


def get_phos_act_inh_matrix(mi: MaudInput):
    """Get two enzyme/phos_enzyme matrix, 1 if enyzme modified, else 0.

    :param mi: a MaudInput object
    """
    enzyme_ix = mi.stan_coords.enzymes
    phos_enz_ix = mi.stan_coords.phos_enzs
    S_phos_act = np.zeros([len(enzyme_ix), len(phos_enz_ix)])
    S_phos_inh = np.zeros([len(enzyme_ix), len(phos_enz_ix)])
    for phos in mi.kinetic_model.phosphorylation:
        i = codify(enzyme_ix)[phos.enzyme_id] - 1
        j = codify(phos_enz_ix)[phos.id] - 1
        if phos.activating:
            S_phos_act[i, j] = 1
        elif phos.inhibiting:
            S_phos_inh[i, j] = 1
    return S_phos_act.T.tolist(), S_phos_inh.T.tolist()


def _get_km_lookup(mi: MaudInput) -> List[List[int]]:
    out = pd.DataFrame(
        0,
        index=mi.stan_coords.mics,
        columns=mi.stan_coords.enzymes + mi.stan_coords.drains,
    )
    for i, (mic_id, enz_id) in enumerate(
        zip(mi.stan_coords.km_mics, mi.stan_coords.km_enzs)
    ):
        out.loc[mic_id, enz_id] = i + 1
    return out.values


def _get_conc_init(mi: MaudInput) -> pd.DataFrame:
    """Get the initial mic concentrations for the ODE solver.

    The initial concentration for a measured mic is the measured value.

    The initial concentration for an unmeasured unbalanced mic in each
    experiment is its prior mean.

    The initial concentration for an unmeasured balanced mics is 0.01

    :param mi: a MaudInput object

    """
    cs = mi.stan_coords
    out = mi.priors.priors_conc_unbalanced.location.reindex(columns=cs.mics).fillna(
        0.01
    )
    for (exp_id, mic_id), value in mi.measurements.yconc["measurement"].items():
        out.loc[exp_id, mic_id] = value
    return out


def get_prior_dict(ps: PriorSet) -> dict:
    """Get a dictionary of priors that be used as a Stan input."""

    def unpack(p: Union[IndPrior1d, IndPrior2d]) -> List[np.ndarray]:
        return [p.location.values.tolist(), p.scale.values.tolist()]

    return {
        "priors_dgf": unpack(ps.priors_dgf),
        "priors_kcat": unpack(ps.priors_kcat),
        "priors_km": unpack(ps.priors_km),
        "priors_ki": unpack(ps.priors_ki),
        "priors_diss_t": unpack(ps.priors_diss_t),
        "priors_diss_r": unpack(ps.priors_diss_r),
        "priors_kcat_phos": unpack(ps.priors_kcat_phos),
        "priors_transfer_constant": unpack(ps.priors_transfer_constant),
        "priors_conc_unbalanced": unpack(ps.priors_conc_unbalanced),
        "priors_conc_enzyme": unpack(ps.priors_conc_enzyme),
        "priors_conc_phos": unpack(ps.priors_conc_phos),
        "priors_drain": unpack(ps.priors_drain),
    }


def get_config_dict(mi: MaudInput) -> dict:
    """Get a dictionary of Stan configuration."""
    config = {
        **{
            "LIKELIHOOD": int(mi.config.likelihood),
            "reject_non_steady": int(mi.config.reject_non_steady),
            "conc_init": _get_conc_init(mi).values,
        },
        **mi.config.ode_config,
    }
    config["max_num_steps"] = int(config["max_num_steps"])
    config["abs_tol_forward"] = [config["abs_tol_forward"]] * len(
        mi.stan_coords.balanced_mics
    )
    config["abs_tol_backward"] = [config["abs_tol_backward"]] * len(
        mi.stan_coords.balanced_mics
    )
    return config


def get_measurements_dict(ms: MeasurementSet, cs: StanCoordSet) -> dict:
    """Get a dictionary of measurements and their indexes."""
    mic_ix = codify(cs.mics)
    rxn_ix = codify(cs.reactions)
    enz_ix = codify(cs.enzymes)
    exp_ix = codify(cs.experiments)
    return {
        "yconc": ms.yconc["measurement"].values.tolist(),
        "sigma_conc": ms.yconc["error_scale"].values.tolist(),
        "experiment_yconc": (
            ms.yconc.reset_index()["experiment_id"].map(exp_ix).values.tolist()
        ),
        "mic_ix_yconc": (
            ms.yconc.reset_index()["target_id"].map(mic_ix).values.tolist()
        ),
        "yflux": ms.yflux["measurement"].values.tolist(),
        "sigma_flux": ms.yflux["error_scale"].values.tolist(),
        "experiment_yflux": (
            ms.yflux.reset_index()["experiment_id"].map(exp_ix).values.tolist()
        ),
        "reaction_yflux": (
            ms.yflux.reset_index()["target_id"].map(rxn_ix).values.tolist()
        ),
        "yenz": ms.yenz["measurement"].values.tolist(),
        "sigma_enz": ms.yenz["error_scale"].values.tolist(),
        "experiment_yenz": (
            ms.yenz.reset_index()["experiment_id"].map(exp_ix).values.tolist()
        ),
        "enzyme_yenz": (ms.yenz.reset_index()["target_id"].map(enz_ix).values.tolist()),
    }


def get_modifier_info(mi: MaudInput) -> Dict[str, pd.Series]:
    """Get a dictionary of modifier info from a MaudInput.

    The output dictionary has the form modifier type -> pandas series indexed
    by edges, whose values are lists of mic ids of the respective modifier.

    For example:

    {
        "competitive_inhibitor": pd.Series({"r2": ["m2"]}, index=["r1", "r2"]),
        "allosteric"_inhibitor: pd.Series({}, index=["r1", "r2"]),
        "allosteric"_activator: pd.Series({"r1": ["m1"]}, index=["r1", "r2"])
    }

    :param mi: A MaudInput object.

    """
    out = {}
    S = get_stoichiometry(mi)
    mods = list(
        itertools.chain(
            *[
                m
                for r in mi.kinetic_model.reactions
                for e in r.enzymes
                for m in e.modifiers.values()
            ]
        )
    )
    for modifier_type in [
        "competitive_inhibitor",
        "allosteric_inhibitor",
        "allosteric_activator",
    ]:
        mod_mics: dict = {edge: [] for edge in S.columns}
        for m in mods:
            if m.modifier_type == modifier_type:
                mod_mics[m.enzyme_id].append(m.mic_id)
        out[modifier_type] = pd.Series(mod_mics).reindex(S.columns)
    return out


def get_phos_info(mi: MaudInput):
    """Get a dictionary of phosphorylation info from a MaudInput.

    The output dictionary has the form phosphorylation type -> pandas series
    indexed by edges, whose values are lists of mic ids of the respective
    modifier.

    For example:

    {
        "activators": pd.Series({"r2": ["m2"]}, index=["r1", "r2"]),
        "inhibitors": pd.Series({"r1": ["m1"]}, index=["r1", "r2"])
    }

    :param mi: A MaudInput object.

    """
    S = get_stoichiometry(mi)
    acts: dict = {edge: [] for edge in S.columns}
    inhs: dict = {edge: [] for edge in S.columns}
    for p in mi.kinetic_model.phosphorylation:
        if p.activating:
            acts[p.enzyme_id].append(p.id)
        elif p.inhibiting:
            inhs[p.enzyme_id].append(p.id)
        else:
            raise ValueError(
                f"Phosphorylation {p.id} is neither inhibiting nor activating."
            )
    return {
        "activators": pd.Series(acts).reindex(S.columns),
        "inhibitors": pd.Series(inhs).reindex(S.columns),
    }


def get_edge_type(edge_id: str, mi: MaudInput):
    """Find the type of an edge, given the id and a MaudInput.

    :param edge_id: string identifying the edge
    """
    reaction_mechanism = [
        rxn.reaction_mechanism
        for rxn in mi.kinetic_model.reactions
        for enz in rxn.enzymes
        if edge_id in enz.id
    ]
    if reaction_mechanism == []:
        reaction_mechanism = [
            rxn.reaction_mechanism
            for rxn in mi.kinetic_model.reactions
            if edge_id in rxn.id
        ]
    reaction_mechanism = reaction_mechanism[0]
    if reaction_mechanism == "reversible_modular_rate_law":
        return 1
    elif reaction_mechanism == "drain":
        return 2
    elif reaction_mechanism == "irreversible_modular_rate_law":
        return 3
    else:
        raise ValueError(f"Edge {reaction_mechanism} is of unknown type.")


def get_input_data(mi: MaudInput) -> dict:
    """Get the input to inference_model.stan from a MaudInput object.

    :param mi: a MaudInput object

    """
    # validate_specified_fluxes(mi)
    sorted_enzymes = sorted(
        [e for r in mi.kinetic_model.reactions for e in r.enzymes],
        key=lambda e: codify(mi.stan_coords.enzymes)[e.id],
    )
    sorted_mics = sorted(
        mi.kinetic_model.mics,
        key=lambda m: codify(mi.stan_coords.mics)[m.id],
    )
    phos_enz_ix = codify(mi.stan_coords.phos_enzs)
    S = get_stoichiometry(mi)
    edge_type = [get_edge_type(eid, mi) for eid in S.columns]
    edge_to_enzyme = (
        pd.Series(
            {e: codify(mi.stan_coords.enzymes)[e] for e in mi.stan_coords.enzymes},
            index=mi.stan_coords.edges,
        )
        .fillna(0)
        .astype(int)
    )
    edge_to_drain = (
        pd.Series(
            {d: codify(mi.stan_coords.drains)[d] for d in mi.stan_coords.drains},
            index=mi.stan_coords.edges,
        )
        .fillna(0)
        .astype(int)
    )
    edge_to_reaction = pd.Series(
        {
            **{
                e.id: codify(mi.stan_coords.reactions)[e.reaction_id]
                for e in sorted_enzymes
            },
            **{d: codify(mi.stan_coords.reactions)[d] for d in mi.stan_coords.drains},
        },
        index=mi.stan_coords.edges,
    ).astype(int)
    phos_info = get_phos_info(mi)
    mod_info = get_modifier_info(mi)
    water_stoichiometry_enzyme = {e.id: e.water_stoichiometry for e in sorted_enzymes}
    water_stoichiometry = pd.Series(water_stoichiometry_enzyme, index=S.columns).fillna(
        0
    )
    mic_to_met = [
        codify(mi.stan_coords.metabolites)[mic.metabolite_id] for mic in sorted_mics
    ]
    knockout_matrix_enzyme = mi.measurements.enz_knockouts.astype(int)
    knockout_matrix_phos = mi.measurements.phos_knockouts.astype(int)
    unbalanced_mic_ix, balanced_mic_ix = (
        [codify(mi.stan_coords.mics)[mic_id] for mic_id in ix]
        for ix in (
            mi.stan_coords.unbalanced_mics,
            mi.stan_coords.balanced_mics,
        )
    )
    allosteric_enzymes = [e for e in sorted_enzymes if e.allosteric]
    prior_dict = get_prior_dict(mi.priors)
    config_dict = get_config_dict(mi)
    measurements_dict = get_measurements_dict(mi.measurements, mi.stan_coords)
    mic_ix = codify(mi.stan_coords.mics)
    return {
        **{
            # sizes
            "N_mic": len(mi.kinetic_model.mics),
            "N_unbalanced": len(mi.stan_coords.unbalanced_mics),
            "N_metabolite": len(mi.stan_coords.metabolites),
            "N_km": len(mi.stan_coords.km_mics),
            "N_reaction": len(mi.stan_coords.reactions),
            "N_enzyme": len(mi.stan_coords.enzymes),
            "N_edge": S.shape[1],
            "N_phosphorylation_enzymes": len(mi.stan_coords.phos_enzs),
            "N_experiment": len(mi.stan_coords.experiments),
            "N_flux_measurement": len(mi.stan_coords.yflux_rxns),
            "N_enzyme_measurement": len(mi.stan_coords.yenz_enzs),
            "N_conc_measurement": len(mi.stan_coords.yconc_mics),
            "N_ki": len(mi.stan_coords.ci_mics),
            "N_ai": len(mi.stan_coords.ai_mics),
            "N_aa": len(mi.stan_coords.aa_mics),
            "N_ae": len(allosteric_enzymes),
            "N_pa": int(phos_info["activators"].apply(lambda l: len(l) > 0).sum()),
            "N_pi": int(phos_info["inhibitors"].apply(lambda l: len(l) > 0).sum()),
            "N_drain": len(mi.stan_coords.drains),
            # indexes
            "unbalanced_mic_ix": unbalanced_mic_ix,
            "balanced_mic_ix": balanced_mic_ix,
            "ix_ci": mod_info["competitive_inhibitor"]
            .explode()
            .map(mic_ix)
            .dropna()
            .astype(int)
            .values,
            "ix_ai": mod_info["allosteric_inhibitor"]
            .explode()
            .map(mic_ix)
            .dropna()
            .astype(int)
            .values,
            "ix_aa": mod_info["allosteric_activator"]
            .explode()
            .map(mic_ix)
            .dropna()
            .astype(int)
            .values,
            "ix_pa": phos_info["activators"]
            .explode()
            .map(phos_enz_ix)
            .dropna()
            .astype(int)
            .values,
            "ix_pi": phos_info["inhibitors"]
            .explode()
            .map(phos_enz_ix)
            .dropna()
            .astype(int)
            .values,
            # network properties
            "S": S.values,
            "edge_type": edge_type,
            "edge_to_enzyme": edge_to_enzyme.values,
            "edge_to_drain": edge_to_drain.values,
            "edge_to_reaction": edge_to_reaction.values,
            "water_stoichiometry": water_stoichiometry.values,
            "mic_to_met": mic_to_met,
            "km_lookup": _get_km_lookup(mi),
            "is_knockout": knockout_matrix_enzyme.values.tolist(),
            "is_phos_knockout": knockout_matrix_phos.values.tolist(),
            "subunits": [e.subunits for e in sorted_enzymes],
            "n_ci": mod_info["competitive_inhibitor"].map(len).values,
            "n_ai": mod_info["allosteric_inhibitor"].map(len).values,
            "n_aa": mod_info["allosteric_activator"].map(len).values,
            "n_pa": phos_info["activators"].map(len).values,
            "n_pi": phos_info["inhibitors"].map(len).values,
        },
        **prior_dict,
        **measurements_dict,
        **config_dict,
    }
