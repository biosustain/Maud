"""Provides function get_stan_variable_set."""

from typing import List, Tuple

from maud.data_model.hardcoding import ID_SEPARATOR
from maud.data_model.kinetic_model import KineticModel, ReactionMechanism
from maud.data_model.measurement_set import MeasurementSet
from maud.data_model.stan_variable_set import (
    ConcEnzyme,
    ConcPhos,
    ConcUnbalanced,
    Dgf,
    DissociationConstant,
    Drain,
    Kcat,
    KcatPhos,
    Ki,
    Km,
    StanVariableSet,
    TransferConstant,
)


AllostericCoords = Tuple[List[str], List[str], List[str], List[str]]
KmCoords = Tuple[List[str], List[str], List[str], List[str]]
KcatCoords = Tuple[List[str], List[str], List[str]]
KiCoords = Tuple[List[str], List[str], List[str], List[str], List[str]]


def get_km_coords(kinetic_model: KineticModel) -> KmCoords:
    km_ids = []
    km_enzs = []
    km_mets = []
    km_cpts = []
    for er in kinetic_model.ers:
        rxn = [r for r in kinetic_model.reactions if r.id == er.reaction_id][0]
        enz = [e for e in kinetic_model.enzymes if e.id == er.enzyme_id][0]
        for mic_id in rxn.stoichiometry.keys():
            km_id = ID_SEPARATOR.join([enz.id, mic_id])
            met_id, cpt_id = mic_id.split(ID_SEPARATOR)
            if km_id not in km_ids:
                km_ids.append(km_id)
                km_enzs.append(enz.id)
                km_mets.append(met_id)
                km_cpts.append(cpt_id)
    return km_ids, km_enzs, km_mets, km_cpts


def get_dc_coords(kinetic_model: KineticModel) -> AllostericCoords:
    ids = []
    enzs = []
    mets = []
    cpts = []
    if kinetic_model.allosteries is not None:
        for a in kinetic_model.allosteries:
            ids.append(a.id)
            enzs.append(a.enzyme_id)
            mets.append(a.metabolite_id)
            cpts.append(a.compartment_id)
    return ids, enzs, mets, cpts


def get_ci_coords(kinetic_model: KineticModel) -> KiCoords:
    ids = []
    enzs = []
    rxns = []
    mets = []
    cpts = []
    if kinetic_model.competitive_inhibitions is not None:
        for ci in kinetic_model.competitive_inhibitions:
            ids.append(ci.id)
            enzs.append(ci.enzyme_id)
            rxns.append(ci.reaction_id)
            mets.append(ci.metabolite_id)
            cpts.append(ci.compartment_id)
    return ids, enzs, rxns, mets, cpts


def get_kcat_coords(kinetic_model: KineticModel) -> KcatCoords:
    ids = []
    enzs = []
    rxns = []
    for er in kinetic_model.ers:
        ids.append(er.id)
        enzs.append(er.enzyme_id)
        rxns.append(er.reaction_id)
    return ids, enzs, rxns


def get_stan_variable_set(kmod: KineticModel, ms: MeasurementSet):
    km_ids, km_enzs, km_mets, km_cpts = get_km_coords(kmod)
    dc_ids, dc_enzs, dc_mets, dc_cpts = get_dc_coords(kmod)
    ci_ids, ci_enzs, ci_rxns, ci_mets, ci_cpts = get_ci_coords(kmod)
    enzyme_ids = [e.id for e in kmod.enzymes]
    kcat_ids, kcat_enzs, kcat_rxns = get_kcat_coords(kmod)
    allosteric_enzyme_ids = (
        list(set(a.enzyme_id for a in kmod.allosteries))
        if kmod.allosteries is not None
        else []
    )
    metabolite_ids = [m.id for m in kmod.metabolites]
    phos_enzs = (
        [p.enzyme_id for p in kmod.phosphorylations]
        if kmod.phosphorylations is not None
        else []
    )
    exp_ids = [e.id for e in ms.experiments]
    drain_ids = [
        d.id for d in kmod.reactions if d.mechanism == ReactionMechanism.DRAIN
    ]
    unbalanced_mics = [m.id for m in kmod.mics if not m.balanced]
    return StanVariableSet(
        dgf=Dgf(ids=[metabolite_ids]),
        km=Km(ids=[km_ids], split_ids=[km_enzs, km_mets, km_cpts]),
        kcat=Kcat(ids=[kcat_ids], split_ids=[kcat_enzs, kcat_rxns]),
        ki=Ki(ids=[ci_ids], split_ids=[ci_enzs, ci_rxns, ci_mets, ci_cpts]),
        dissociation_constant=DissociationConstant(ids=[dc_ids], split_ids=[dc_enzs, dc_mets, dc_cpts]),
        transfer_constant=TransferConstant(ids=[allosteric_enzyme_ids]),
        kcat_phos=KcatPhos(ids=[phos_enzs]),
        drain=Drain(ids=[exp_ids, drain_ids]),
        conc_enzyme=ConcEnzyme(ids=[exp_ids, enzyme_ids]),
        conc_unbalanced=ConcUnbalanced(ids=[exp_ids, unbalanced_mics]),
        conc_phos=ConcPhos(ids=[exp_ids, phos_enzs]),
    )
