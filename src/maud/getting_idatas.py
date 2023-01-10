"""Functions for creating InferenceData objects from Maud outputs."""
from typing import List

import arviz as az

from maud.data_model.experiment import MeasurementType
from maud.data_model.hardcoding import ID_SEPARATOR
from maud.data_model.maud_input import MaudInput


def get_idata(csvs: List[str], mi: MaudInput, mode: str) -> az.InferenceData:
    """Get an arviz InferenceData object from Maud csvs."""

    experiments = (
        [e for e in mi.experiments if e.is_train]
        if mode == "train"
        else [e for e in mi.experiments if e.is_test]
    )
    yconc_coords, yflux_coords, yenz_coords = (
        [
            f"{e.id}{ID_SEPARATOR}{m.target_id}"
            for e in experiments
            for m in e.measurements
            if m.target_type == t
        ]
        for t in [
            MeasurementType.MIC,
            MeasurementType.FLUX,
            MeasurementType.ENZYME,
        ]
    )
    coords = {
        "enzymes": [e.id for e in mi.kinetic_model.enzymes],
        "experiments": [e.id for e in experiments],
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
            pme.id for pme in mi.kinetic_model.phosphorylation_modifying_enzymes
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
        "kms": mi.parameters.km.ids[0],
        "kis": mi.parameters.ki.ids[0],
        "dissociation_constants": (mi.parameters.dissociation_constant.ids[0]),
        "yconcs": yconc_coords,
        "yfluxs": yflux_coords,
        "yenz": yenz_coords,
    }
    dims = {
        f"flux_{mode}": ["experiments", "reactions"],
        f"conc_{mode}": ["experiments", "mics"],
        f"log_conc_enzyme_{mode}_z": ["experiments", "enzymes"],
        f"conc_enzyme_{mode}": ["experiments", "enzymes"],
        f"conc_unbalanced_{mode}": ["experiments", "unbalanced_mics"],
        f"conc_pme_{mode}": [
            "experiments",
            "phosphorylation_modifying_enzymes",
        ],
        f"drain_{mode}": ["experiments", "drains"],
        f"psi_{mode}": ["experiments"],
        f"saturation_{mode}": ["experiments", "edges"],
        f"allostery_{mode}": ["experiments", "edges"],
        f"phosphorylation_{mode}": ["experiments", "edges"],
        f"reversibility_{mode}": ["experiments", "edges"],
        "dissociation_constant": ["allosteries"],
        "log_transfer_constant_z": ["allosteric_enzymes"],
        "transfer_constant": ["allosteric_enzymes"],
        "dgf": ["metabolites"],
        "dgr_train": ["experiments", "edges"],
        "keq": ["experiments", "edges"],
        "kcat": ["enzymes"],
        "kcat_pme": ["phosphorylation_modifying_enzymes"],
        "km": ["kms"],
        "ki": ["kis"],
        f"yrep_conc_{mode}": ["yconcs"],
        f"yrep_flux_{mode}": ["yfluxs"],
        f"llik_conc_{mode}": ["yconcs"],
        f"llik_flux_{mode}": ["yfluxs"],
    }
    # for zvar in [
    #     "kcat",
    #     "km",
    #     "kcat_pme",
    #     "ki",
    #     "dissociation_constant",
    #     "transfer_constant",
    #     "psi_train",
    #     f"drain_{mode}",
    #     f"conc_unbalanced_{mode}",
    #     f"conc_pme_{mode}"
    # ]:
    #     dims[f"log_{zvar}_z"] = dims[zvar]
    return az.from_cmdstan(
        posterior=csvs,
        coords=coords,
        log_likelihood=[f"llik_conc_{mode}", f"llik_flux_{mode}"],
        posterior_predictive=[f"yrep_conc_{mode}", f"yrep_flux_{mode}"],
        dims=dims,
        save_warmup=True,
    )
