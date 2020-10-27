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

import os
from copy import deepcopy
from itertools import product
from typing import Dict

import cmdstanpy
import numpy as np
import pandas as pd

from maud import io
from maud.data_model import KineticModel, MaudInput


HERE = os.path.dirname(os.path.abspath(__file__))
INCLUDE_PATH = ""
DEFAULT_PRIOR_LOC_UNBALANCED = 0.1
DEFAULT_PRIOR_SCALE_UNBALANCED = 2
DEFAULT_PRIOR_LOC_ENZYME = 0.1
DEFAULT_PRIOR_SCALE_ENZYME = 2
STAN_PROGRAM_RELATIVE_PATH = "inference_model.stan"


def get_full_stoichiometry(
    kinetic_model: KineticModel,
    enzyme_codes: Dict[str, int],
    metabolite_codes: Dict[str, int],
    reaction_codes: Dict[str, int],
    drain_codes,
):
    """Get full stoichiometric matrix for each isoenzyme.

    :param kinetic_model: A Kinetic Model object
    :param enzyme_codes: the codified enzyme codes
    :param metabolite_codes: the codified metabolite codes
    """
    S_enz = pd.DataFrame(0, index=enzyme_codes, columns=metabolite_codes)
    if drain_codes is not None:
        S_drain = pd.DataFrame(0, index=drain_codes, columns=metabolite_codes)
        S_complete = pd.DataFrame(
            0, index={**enzyme_codes, **drain_codes}, columns=metabolite_codes
        )
    else:
        S_complete = pd.DataFrame(0, index={**enzyme_codes}, columns=metabolite_codes)
    S_enz_to_flux_map = pd.DataFrame(0, index=reaction_codes, columns=enzyme_codes)

    for rxn_id, rxn in kinetic_model.reactions.items():
        for enz_id, _ in rxn.enzymes.items():
            for met, stoic in rxn.stoichiometry.items():
                S_enz_to_flux_map.loc[rxn_id, enz_id] = 1
                S_enz.loc[enz_id, met] = stoic
                S_complete.loc[enz_id, met] = stoic

    if drain_codes is not None:
        for drain_id, drain in kinetic_model.drains.items():
            for met, stoic in drain.stoichiometry.items():
                S_drain.loc[drain_id, met] = stoic
                S_complete.loc[drain_id, met] = stoic

    if drain_codes is not None:
        return S_enz, S_enz_to_flux_map, S_complete, S_drain
    else:
        return S_enz, S_enz_to_flux_map, S_complete


def get_knockout_matrix(mi: MaudInput):
    """Get binary experiment, enzyme matrix, 1 if enyzme knocked out, 0 if not.

    :param mi: a MaudInput object
    """
    experiment_codes = mi.stan_codes["experiment"]
    enzyme_codes = mi.stan_codes["enzyme"]
    enzyme_knockout_matrix = pd.DataFrame(
        0, index=np.arange(len(experiment_codes)), columns=np.arange(len(enzyme_codes))
    )
    for exp_name, exp in mi.experiments.items():
        if exp.knockouts:
            for enz in exp.knockouts:
                enzyme_knockout_matrix.loc[
                    experiment_codes[exp_name] - 1, enzyme_codes[enz] - 1
                ] = 1
    return enzyme_knockout_matrix


def get_drain_priors(mi: MaudInput):
    """Returns an experiment x drain matrix of location and scale priors for drains."""
    drain_codes = mi.stan_codes["drain"]
    enzyme_codes = mi.stan_codes["enzyme"]
    experiment_codes = mi.stan_codes["experiment"]

    drain_loc_prior = pd.DataFrame(index=experiment_codes, columns=drain_codes)
    drain_scale_prior = pd.DataFrame(index=experiment_codes, columns=drain_codes)

    for p in mi.priors["drains"]:
        drain_loc_prior.loc[p.experiment_id, p.drain_id] = p.location
        drain_scale_prior.loc[p.experiment_id, p.drain_id] = p.scale
    return drain_loc_prior, drain_scale_prior


def get_km_lookup(km_priors, mic_codes, enzyme_codes):
    """Get a mics x enzymes matrix of km indexes. 0 means null."""
    out = pd.DataFrame(0, index=mic_codes.values(), columns=enzyme_codes.values())
    for i, row in km_priors.iterrows():
        mic_code = mic_codes[row["mic_id"]]
        enzyme_code = enzyme_codes[row["enzyme_id"]]
        out.loc[mic_code, enzyme_code] = i + 1
    return out


def sample(
    data_path: str,
    rel_tol: float,
    abs_tol: float,
    max_num_steps: int,
    likelihood: int,
    n_samples: int,
    n_warmup: int,
    n_chains: int,
    n_cores: int,
    timepoint: float,
    output_dir: str,
    threads_per_chain: int,
    save_warmup: bool,
) -> cmdstanpy.CmdStanMCMC:
    """Sample from a posterior distribution.

    :param data_path: A path to a toml file containing input data
    :param rel_tol: Sets ODE solver's rel_tol control parameter
    :param abs_tol: Sets ODE solver's abs_tol control parameter
    :param max_num_steps: Sets ODE solver's max_num_steps control parameter
    :param likelihood: Set to zero if you don't want the model to include information
    from experimental data.
    :param n_samples: Number of post-warmup samples
    :param n_warmup: Number of warmup samples
    :param n_chains: Number of MCMC chains to run
    :param n_cores: Number of cores to try and use
    :param time_step: Amount of time for the ode solver to simulate in order to compare
    initial state with evolved state
    :param: output_dir: Directory to save output
    :param: threads_per_chain: Number of threads per chain (default is 1)
    :param: save_warmup: whether to save warmup draws (default is True)
    """
    if threads_per_chain != 1:
        os.environ["STAN_NUM_THREADS"] = str(threads_per_chain)
    model_name = os.path.splitext(os.path.basename(data_path))[0]
    input_filepath = os.path.join(output_dir, f"input_data_{model_name}.json")
    init_filepath = os.path.join(output_dir, f"initial_values_{model_name}.json")
    mi = io.load_maud_input_from_toml(data_path)
    input_data = get_input_data(
        mi, abs_tol, rel_tol, max_num_steps, likelihood, timepoint
    )
    inits = get_initial_conditions(input_data, mi)
    cmdstanpy.utils.jsondump(input_filepath, input_data)
    cmdstanpy.utils.jsondump(init_filepath, inits)
    stan_program_filepath = os.path.join(HERE, STAN_PROGRAM_RELATIVE_PATH)
    include_path = os.path.join(HERE, INCLUDE_PATH)
    cpp_options = {"STAN_THREADS": True} if threads_per_chain != 1 else None
    model = cmdstanpy.CmdStanModel(
        stan_file=stan_program_filepath,
        stanc_options={"include_paths": [include_path]},
        cpp_options=cpp_options,
    )
    return model.sample(
        data=input_filepath,
        chains=n_chains,
        parallel_chains=n_cores,
        iter_sampling=n_samples,
        output_dir=output_dir,
        iter_warmup=n_warmup,
        max_treedepth=11,
        save_warmup=save_warmup,
        inits=init_filepath,
        show_progress=True,
        step_size=0.025,
        adapt_delta=0.99,
    )


def simulate_once(
    data_path: str,
    rel_tol: float,
    abs_tol: float,
    max_num_steps: float,
    timepoint: float,
    output_dir: str,
) -> cmdstanpy.CmdStanMCMC:
    """Sample from a posterior distribution.

    :param data_path: A path to a toml file containing input data
    :param rel_tol: Sets ODE solver's rel_tol control parameter
    :param abs_tol: Sets ODE solver's abs_tol control parameter
    :param max_num_steps: Sets ODE solver's max_num_steps control parameter
    from experimental data.
    initial state with evolved state
    :param: output_dir: Directory to save output
    """
    model_name = os.path.splitext(os.path.basename(data_path))[0]
    input_filepath = os.path.join(output_dir, f"input_data_{model_name}.json")
    init_filepath = os.path.join(output_dir, f"initial_values_{model_name}.json")
    mi = io.load_maud_input_from_toml(data_path)
    input_data = get_input_data(mi, abs_tol, rel_tol, max_num_steps, 1, timepoint)
    inits = get_initial_conditions(input_data, mi)
    cmdstanpy.utils.jsondump(input_filepath, input_data)
    cmdstanpy.utils.jsondump(init_filepath, inits)
    stan_program_filepath = os.path.join(HERE, STAN_PROGRAM_RELATIVE_PATH)
    include_path = os.path.join(HERE, INCLUDE_PATH)
    model = cmdstanpy.CmdStanModel(
        stan_file=stan_program_filepath, stanc_options={"include_paths": [include_path]}
    )
    return model.sample(
        data=input_filepath,
        chains=1,
        iter_sampling=1,
        output_dir=output_dir,
        inits=init_filepath,
        show_progress=True,
        fixed_param=True,
    )


def get_input_data(
    mi: MaudInput,
    abs_tol: float,
    rel_tol: float,
    max_num_steps: int,
    likelihood: int,
    timepoint: float,
) -> dict:
    """Put a MaudInput and some config numbers into a Stan-friendly dictionary.

    :param mi: a MaudInput object
    :param abs_tol: Sets ODE solver's abs_tol control parameter
    :param rel_tol: Sets ODE solver's rel_tol control parameter
    :param max_num_steps: Sets ODE solver's max_num_steps control parameter
    :param likelihood: Set to zero if you don't want the model to include information
    """
    metabolites = mi.kinetic_model.metabolites
    mics = mi.kinetic_model.mics
    reactions = mi.kinetic_model.reactions
    reaction_codes = mi.stan_codes["reaction"]
    drain_codes = mi.stan_codes["drain"]
    enzyme_codes = mi.stan_codes["enzyme"]
    experiment_codes = mi.stan_codes["experiment"]
    met_codes = mi.stan_codes["metabolite"]
    mic_codes = mi.stan_codes["metabolite_in_compartment"]
    balanced_mic_codes = mi.stan_codes["balanced_mic"]
    unbalanced_mic_codes = mi.stan_codes["unbalanced_mic"]
    mic_to_met = {
        mic_codes[mic.id]: met_codes[mic.metabolite_id] for mic in mics.values()
    }
    enzymes = {k: v for r in reactions.values() for k, v in r.enzymes.items()}
    water_stoichiometry = [
        r.water_stoichiometry for r in reactions.values() for e in r.enzymes.items()
    ]
    if drain_codes is not None:
        (
            enzyme_stoic,
            stoic_map_to_flux,
            full_stoic,
            drain_stoic,
        ) = get_full_stoichiometry(
            mi.kinetic_model, enzyme_codes, mic_codes, reaction_codes, drain_codes
        )
    else:
        enzyme_stoic, stoic_map_to_flux, full_stoic = get_full_stoichiometry(
            mi.kinetic_model, enzyme_codes, mic_codes, reaction_codes, drain_codes
        )
        drain_stoic = pd.DataFrame()
    subunits = pd.DataFrame(
        {
            "enzyme_id": [e.id for e in enzymes.values()],
            "subunits": [e.subunits for e in enzymes.values()],
            "index": [enzyme_codes[e.id] for e in enzymes.values()],
        }
    )

    subunits = subunits.sort_values(by="index").reset_index(drop=True)

    # priors
    km_priors = pd.DataFrame(
        {
            "location": [p.location for p in mi.priors["kms"]],
            "scale": [p.scale for p in mi.priors["kms"]],
            "enzyme_id": [p.enzyme_id for p in mi.priors["kms"]],
            "mic_id": [p.mic_id for p in mi.priors["kms"]],
        }
    )
    kcat_priors, transfer_constant_priors = (
        pd.DataFrame(
            {
                "location": [p.location for p in mi.priors[prior_type]],
                "scale": [p.scale for p in mi.priors[prior_type]],
                "enzyme_id": [p.enzyme_id for p in mi.priors[prior_type]],
            }
        )
        for prior_type in ["kcats", "transfer_constants"]
    )
    formation_energy_priors = pd.DataFrame(
        {
            "location": [p.location for p in mi.priors["formation_energies"]],
            "scale": [p.scale for p in mi.priors["formation_energies"]],
            "metabolite_id": [p.metabolite_id for p in mi.priors["formation_energies"]],
        }
    )
    formation_energy_priors["stan_code"] = formation_energy_priors["metabolite_id"].map(
        met_codes
    )
    formation_energy_priors.sort_values("stan_code", inplace=True)
    if drain_codes is not None:
        prior_loc_drain, prior_scale_drain = get_drain_priors(mi)
    else:
        prior_loc_drain = pd.DataFrame()
        prior_scale_drain = pd.DataFrame()
    ki_priors, diss_t_priors, diss_r_priors = (
        pd.DataFrame(
            {
                "location": [p.location for p in mi.priors[prior_type]],
                "scale": [p.scale for p in mi.priors[prior_type]],
                "mic_id": [p.mic_id for p in mi.priors[prior_type]],
                "enzyme_id": [p.enzyme_id for p in mi.priors[prior_type]],
            }
        )
        if prior_type in mi.priors.keys()
        else pd.DataFrame()
        for prior_type in [
            "inhibition_constants",
            "tense_dissociation_constants",
            "relaxed_dissociation_constants",
        ]
    )
    n_modifier = {}
    for df, param_type in zip(
        [ki_priors, diss_t_priors, diss_r_priors], ["ki", "diss_t", "diss_r"]
    ):
        df["mic_code"] = df["mic_id"].map(mic_codes)
        df["enzyme_code"] = df["enzyme_id"].map(enzyme_codes)
        df.sort_values("enzyme_code", inplace=True)
        n_modifier[param_type] = (
            df.groupby("enzyme_code")
            .size()
            .reindex(enzyme_codes.values())
            .fillna(0)
            .astype(int)
            .tolist()
        )
    unb_shape = len(mi.experiments), len(unbalanced_mic_codes)
    prior_loc_unb = np.full(unb_shape, DEFAULT_PRIOR_LOC_UNBALANCED)
    prior_scale_unb = np.full(unb_shape, DEFAULT_PRIOR_SCALE_UNBALANCED)
    if "mic_concentrations" in mi.priors.keys():
        for p in mi.priors["mic_concentrations"]:
            ix = [experiment_codes[mic_codes[p.mic_id - 1], p.experiment_id] - 1]
            prior_loc_unb[ix] = p.location
            prior_scale_unb[ix] = p.scale
    enzyme_shape = len(mi.experiments), len(enzymes)
    prior_loc_enzyme = np.full(enzyme_shape, DEFAULT_PRIOR_LOC_ENZYME)
    prior_scale_enzyme = np.full(enzyme_shape, DEFAULT_PRIOR_SCALE_ENZYME)

    # measurements
    mic_measurements, reaction_measurements, enzyme_measurements = (
        pd.DataFrame(
            [
                [exp.id, meas.target_id, meas.value, meas.uncertainty]
                for exp in mi.experiments.values()
                for meas in exp.measurements[measurement_type].values()
            ],
            columns=["experiment_id", "target_id", "value", "uncertainty"],
        )
        for measurement_type in ["metabolite", "reaction", "enzyme"]
    )

    conc_init = pd.DataFrame(
        0.01, index=experiment_codes.values(), columns=mic_codes.values()
    )
    for _i, row in mic_measurements.iterrows():
        if row["target_id"] in balanced_mic_codes.keys():
            row_ix = experiment_codes[row["experiment_id"]]
            column_ix = mic_codes[row["target_id"]]
            conc_init.loc[row_ix, column_ix] = row["value"]
    knockout_matrix = get_knockout_matrix(mi=mi)
    km_lookup = get_km_lookup(km_priors, mic_codes, enzyme_codes)
    return {
        "N_mic": len(mics),
        "N_unbalanced": len(unbalanced_mic_codes),
        "N_metabolite": len(metabolites),
        "N_km": len(km_priors),
        "N_reaction": len(reactions),
        "N_enzyme": len(enzymes),
        "N_experiment": len(mi.experiments),
        "N_flux_measurement": len(reaction_measurements),
        "N_enzyme_measurement": len(enzyme_measurements),
        "N_conc_measurement": len(mic_measurements),
        "N_competitive_inhibitor": len(ki_priors),
        "N_allosteric_inhibitor": len(diss_t_priors),
        "N_allosteric_activator": len(diss_r_priors),
        "N_allosteric_enzyme": len(transfer_constant_priors),
        "N_drain": len(mi.kinetic_model.drains),
        "unbalanced_mic_ix": list(unbalanced_mic_codes.values()),
        "balanced_mic_ix": list(balanced_mic_codes.values()),
        "experiment_yconc": (
            mic_measurements["experiment_id"].map(experiment_codes).values
        ),
        "mic_ix_yconc": mic_measurements["target_id"].map(mic_codes).values,
        "yconc": mic_measurements["value"].values,
        "sigma_conc": mic_measurements["uncertainty"].values,
        "experiment_yflux": (
            reaction_measurements["experiment_id"].map(experiment_codes).values
        ),
        "reaction_yflux": (
            reaction_measurements["target_id"].map(reaction_codes).values
        ),
        "yflux": reaction_measurements["value"].values,
        "sigma_flux": reaction_measurements["uncertainty"].values,
        "experiment_yenz": (
            enzyme_measurements["experiment_id"].map(experiment_codes).values
        ),
        "enzyme_yenz": (enzyme_measurements["target_id"].map(enzyme_codes).values),
        "yenz": enzyme_measurements["value"].values,
        "sigma_enz": enzyme_measurements["uncertainty"].values,
        "prior_loc_formation_energy": formation_energy_priors["location"].values,
        "prior_scale_formation_energy": formation_energy_priors["scale"].values,
        "prior_loc_kcat": kcat_priors["location"].values,
        "prior_scale_kcat": kcat_priors["scale"].values,
        "prior_loc_km": km_priors["location"].values,
        "prior_scale_km": km_priors["scale"].values,
        "prior_loc_ki": ki_priors["location"].values,
        "prior_scale_ki": ki_priors["scale"].values,
        "prior_loc_diss_t": diss_t_priors["location"].values,
        "prior_scale_diss_t": diss_t_priors["scale"].values,
        "prior_loc_diss_r": diss_r_priors["location"].values,
        "prior_scale_diss_r": diss_r_priors["scale"].values,
        "prior_loc_tc": transfer_constant_priors["location"].values,
        "prior_scale_tc": transfer_constant_priors["scale"].values,
        "prior_loc_unbalanced": prior_loc_unb,
        "prior_scale_unbalanced": prior_scale_unb,
        "prior_loc_enzyme": prior_loc_enzyme,
        "prior_scale_enzyme": prior_scale_enzyme,
        "prior_loc_drain": prior_loc_drain.values,
        "prior_scale_drain": prior_scale_drain.values,
        "S_enz": enzyme_stoic.T.values,
        "S_drain": drain_stoic.T.values,
        "S_full": full_stoic.T.values,
        "water_stoichiometry": water_stoichiometry,
        "mic_to_met": list(mic_to_met.values()),
        "S_to_flux_map": stoic_map_to_flux.values,
        "is_knockout": knockout_matrix.values,
        "km_lookup": km_lookup.values,
        "n_ci": n_modifier["ki"],
        "n_ai": n_modifier["diss_t"],
        "n_aa": n_modifier["diss_r"],
        "ci_ix": ki_priors["mic_code"].values,
        "ai_ix": diss_t_priors["mic_code"].values,
        "aa_ix": diss_r_priors["mic_code"].values,
        "subunits": subunits["subunits"].values,
        "conc_init": conc_init.values,
        "rel_tol": rel_tol,
        "abs_tol": abs_tol,
        "max_num_steps": max_num_steps,
        "LIKELIHOOD": likelihood,
        "timepoint": timepoint,
    }


def get_init_enzyme(mi):
    """Get initial enyme concentrations.

    Some internal logic is required to ensure that enzyme concentrations are
    always initialised implied vmax is always at least as high as the measured
    flux.

    :param mi: a MaudInput object

    """

    def get_naive_init(mi):
        enz_codes = mi.stan_codes["enzyme"]
        exp_codes = mi.stan_codes["experiment"]
        exp_ids, enz_ids = zip(*product(exp_codes.keys(), enz_codes.keys()))
        out = pd.DataFrame(
            {
                "experiment_id": exp_ids,
                "enzyme_id": enz_ids,
                "init": DEFAULT_PRIOR_LOC_ENZYME,
            }
        ).set_index(["experiment_id", "enzyme_id"])
        for (exp_id, enz_id), _ in out.iterrows():
            measurements = mi.experiments[exp_id].measurements
            if "enzyme" in measurements.keys():
                if enz_id in measurements["enzyme"].keys():
                    out.loc[(exp_id, enz_id), "init"] = measurements["enzyme"][
                        enz_id
                    ].value
        return out.reset_index()

    def get_greatest_measured_flux(exp_id, enz_id, mi):
        measurements = mi.experiments[exp_id].measurements["reaction"]
        rxn_ids = [
            rid
            for rid, r in mi.kinetic_model.reactions.items()
            if enz_id in r.enzymes.keys() and rid in measurements.keys()
        ]
        return (
            max([measurements[rid].value for rid in rxn_ids])
            if len(rxn_ids) > 0
            else np.nan
        )

    def get_vmax_component(enzyme_id, init, mi):
        kcat = [p.location for p in mi.priors["kcats"] if p.enzyme_id == enzyme_id][0]
        return init * kcat

    def get_lowest_vmax(exp_id, enz_id, inits, mi):
        rxns = [
            r for r in mi.kinetic_model.reactions.values() if enz_id in r.enzymes.keys()
        ]
        vmaxes = []
        for rxn in rxns:
            vmax_components = []
            for enz_id in rxn.enzymes.keys():
                init = inits.set_index(["experiment_id", "enzyme_id"]).loc[
                    (exp_id, enz_id), "init"
                ]
                vmax_components.append(get_vmax_component(enz_id, init, mi))
            vmaxes.append(sum(vmax_components))
        return min(vmaxes)

    e = get_naive_init(mi)
    e["experiment_id_stan"] = e["experiment_id"].map(mi.stan_codes["experiment"])
    e["enzyme_id_stan"] = e["enzyme_id"].map(mi.stan_codes["enzyme"])
    e["measured_flux"] = e.apply(
        lambda row: get_greatest_measured_flux(
            row["experiment_id"], row["enzyme_id"], mi
        ),
        axis=1,
    )
    e["vmax"] = e.apply(
        lambda row: get_lowest_vmax(row["experiment_id"], row["enzyme_id"], e, mi),
        axis=1,
    )
    # replace init with init * 1.5 * flux/vmax if vmax is too low
    e["init_out"] = np.where(
        e["vmax"] < e["measured_flux"],
        e["init"] * 1.5 * e["measured_flux"] / e["vmax"],
        e["init"],
    )
    return (
        e.set_index(["experiment_id_stan", "enzyme_id_stan"])["init_out"]
        .sort_index()
        .unstack()
    )


def get_initial_conditions(input_data, mi):
    """Specify parameters' initial conditions."""
    init_conc_unb = deepcopy(input_data["prior_loc_unbalanced"])
    init_unbalanced = pd.DataFrame(
        init_conc_unb,
        index=range(1, input_data["N_experiment"] + 1),
        columns=input_data["unbalanced_mic_ix"],
    )
    for exp_ix, mic_ix, measurement in zip(
        input_data["experiment_yconc"], input_data["mic_ix_yconc"], input_data["yconc"]
    ):
        if mic_ix in input_data["unbalanced_mic_ix"]:
            init_unbalanced.loc[exp_ix, mic_ix] = measurement
    init_enzyme = get_init_enzyme(mi)
    return {
        "km": input_data["prior_loc_km"],
        "ki": input_data["prior_loc_ki"],
        "diss_t": input_data["prior_loc_diss_t"],
        "diss_r": input_data["prior_loc_diss_r"],
        "transfer_constant": input_data["prior_loc_tc"],
        "kcat": input_data["prior_loc_kcat"],
        "conc_unbalanced": init_unbalanced.values,
        "enzyme": init_enzyme.values,
        "formation_energy": input_data["prior_loc_formation_energy"],
        "drain": input_data["prior_loc_drain"],
    }
