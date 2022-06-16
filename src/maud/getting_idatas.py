"""Functions for creating InferenceData objects from Maud outputs."""

from typing import List

import arviz as az
from arviz.utils import Numba

from maud.data_model.hardcoding import ID_SEPARATOR
from maud.data_model.maud_input import MaudInput
from maud.utils import join_str_cols


def get_idata(csvs: List[str], mi: MaudInput, mode: str) -> az.InferenceData:
    """Get an arviz InferenceData object from Maud csvs."""

    Numba.disable_numba()
    experiments = (
        [e.id for e in mi.measurements.experiments if e.is_train]
        if mode == "train"
        else [e.id for e in mi.measurements.experiments if e.is_test]
    )
    coords = {
        "enzymes": [e.id for e in mi.kinetic_model.enzymes],
        "experiments": experiments,
        "reactions": [r.id for r in mi.kinetic_model.reactions],
        "drains": [r.id for r in mi.kinetic_model.drains],
        "metabolites": [m.id for m in mi.kinetic_model.metabolites],
        "mics": [m.id for m in mi.kinetic_model.mics],
        "edges": [e.id for e in mi.kinetic_model.edges],
        "unbalanced_mics": [
            m.id for m in mi.kinetic_model.mics if not m.balanced
        ],
        "phosphorylations": [p.id for p in mi.kinetic_model.phosphorylations]
        if mi.kinetic_model.phosphorylations is not None
        else [],
        "phosphorylation_modifying_enzymes": [
            pme.id
            for pme in mi.kinetic_model.phosphorylation_modifying_enzymes
        ]
        if mi.kinetic_model.phosphorylation_modifying_enzymes is not None
        else [],
        "allosteries": [p.id for p in mi.kinetic_model.allosteries]
        if mi.kinetic_model.allosteries is not None
        else [],
        "allosteric_enzymes": [
            e.id for e in mi.kinetic_model.allosteric_enzymes
        ]
        if mi.kinetic_model.allosteric_enzymes is not None
        else [],
        "competitive_inhibitions": [
            p.id for p in mi.kinetic_model.competitive_inhibitions
        ]
        if mi.kinetic_model.competitive_inhibitions is not None
        else [],
        "kms": mi.stan_variable_set.km.ids[0],
        "kis": mi.stan_variable_set.ki.ids[0],
        "dissociation_constants": (
            mi.stan_variable_set.dissociation_constant.ids[0]
        ),
        "yconcs": join_str_cols(
            mi.measurements.yconc[["experiment_id", "target_id"]].loc[
                lambda df: df["experiment_id"].isin(experiments)
            ],
            sep=ID_SEPARATOR,
        ).to_list(),
        "yfluxs": join_str_cols(
            mi.measurements.yflux[["experiment_id", "target_id"]].loc[
                lambda df: df["experiment_id"].isin(experiments)
            ],
            sep=ID_SEPARATOR,
        ).to_list(),
        "yenz": join_str_cols(
            mi.measurements.yenz[["experiment_id", "target_id"]].loc[
                lambda df: df["experiment_id"].isin(experiments)
            ],
            sep=ID_SEPARATOR,
        ).to_list(),
    }
    return az.from_cmdstan(
        csvs,
        coords=coords,
        dims={
            "flux": ["experiments", "reactions"],
            "conc": ["experiments", "mics"],
            "conc_enzyme": ["experiments", "enzymes"],
            "conc_unbalanced": ["experiments", "unbalanced_mics"],
            "conc_pme": ["experiments", "phosphorylation_modifying_enzymes"],
            "drain": ["experiments", "drains"],
            "dissociation_constant": ["allosteries"],
            "transfer_constant": ["allosteric_enzymes"],
            "dgf": ["metabolites"],
            "dgrs": ["experiments", "edges"],
            "keq": ["experiments", "edges"],
            "kcat": ["enzymes"],
            "kcat_pme": ["phosphorylation_modifying_enzymes"],
            "km": ["kms"],
            "ki": ["kis"],
            "psi": ["experiments"],
            "yconc_sim": ["yconcs"],
            "yflux_sim": ["yfluxs"],
            "yenz_sim": ["yenzs"],
            "log_lik_conc": ["yconcs"],
            "log_lik_flux": ["yfluxs"],
            "log_lik_enz": ["yenzs"],
            "saturation": ["experiments", "edges"],
            "allostery": ["experiments", "edges"],
            "phosphorylation": ["experiments", "edges"],
            "reversibility": ["experiments", "edges"],
        },
        save_warmup=True,
    )
