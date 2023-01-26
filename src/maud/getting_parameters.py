"""Provides function get_maud_parameters."""

from typing import List, Tuple

from maud.data_model.experiment import Experiment, MeasurementType
from maud.data_model.hardcoding import ID_SEPARATOR
from maud.data_model.kinetic_model import KineticModel, ReactionMechanism
from maud.data_model.maud_init import InitInput
from maud.data_model.maud_parameter import (
    ConcEnzymeTest,
    ConcEnzymeTrain,
    ConcPmeTest,
    ConcPmeTrain,
    ConcUnbalancedTest,
    ConcUnbalancedTrain,
    Dgf,
    DissociationConstant,
    DrainTest,
    DrainTrain,
    Kcat,
    KcatPme,
    Ki,
    Km,
    ParameterSet,
    PsiTest,
    PsiTrain,
    TransferConstant,
)
from maud.data_model.prior_input import PriorInput


AllostericCoords = Tuple[List[str], List[str], List[str], List[str], List[str]]
KmCoords = Tuple[List[str], List[str], List[str], List[str]]
KcatCoords = Tuple[List[str], List[str], List[str]]
KiCoords = Tuple[List[str], List[str], List[str], List[str], List[str]]


def get_km_coords(kinetic_model: KineticModel) -> KmCoords:
    """Get ids and split ids for a model's Kms."""
    km_ids = []
    km_enzs = []
    km_mets = []
    km_cpts = []
    for er in kinetic_model.ers:
        rxn = [r for r in kinetic_model.reactions if r.id == er.reaction_id][0]
        enz = [e for e in kinetic_model.enzymes if e.id == er.enzyme_id][0]
        mic_ids = (
            list(rxn.stoichiometry.keys())
            if rxn.mechanism != ReactionMechanism.irreversible_michaelis_menten
            else [k for k, v in rxn.stoichiometry.items() if v < 0]
        )
        for mic_id in mic_ids:
            km_id = ID_SEPARATOR.join([enz.id, mic_id])
            met_id, cpt_id = mic_id.split(ID_SEPARATOR)
            if km_id not in km_ids:
                km_ids.append(km_id)
                km_enzs.append(enz.id)
                km_mets.append(met_id)
                km_cpts.append(cpt_id)
    return km_ids, km_enzs, km_mets, km_cpts


def get_dc_coords(kinetic_model: KineticModel) -> AllostericCoords:
    """Get ids and split ids for a model's dissociation constants."""
    ids = []
    enzs = []
    mets = []
    cpts = []
    mts = []
    if kinetic_model.allosteries is not None:
        for a in kinetic_model.allosteries:
            ids.append(a.id)
            enzs.append(a.enzyme_id)
            mets.append(a.metabolite_id)
            cpts.append(a.compartment_id)
            mts.append(a.modification_type.name)
    return ids, enzs, mets, cpts, mts


def get_ci_coords(kinetic_model: KineticModel) -> KiCoords:
    """Get ids and split ids for a model's competitive inhibition constants."""
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
    """Get ids and split ids for a model's Kcats."""
    ids = []
    enzs = []
    rxns = []
    for er in kinetic_model.ers:
        ids.append(er.id)
        enzs.append(er.enzyme_id)
        rxns.append(er.reaction_id)
    return ids, enzs, rxns


def get_maud_parameters(
    kmod: KineticModel,
    experiments: List[Experiment],
    pi: PriorInput,
    ii: InitInput,
):
    """Get a ParameterSet object from a KineticModel and a MeasurementSet."""
    km_ids, km_enzs, km_mets, km_cpts = get_km_coords(kmod)
    dc_ids, dc_enzs, dc_mets, dc_cpts, dc_mts = get_dc_coords(kmod)
    ci_ids, ci_enzs, ci_rxns, ci_mets, ci_cpts = get_ci_coords(kmod)
    enzyme_ids = [e.id for e in kmod.enzymes]
    kcat_ids, kcat_enzs, kcat_rxns = get_kcat_coords(kmod)
    allosteric_enzyme_ids = (
        [e.id for e in kmod.allosteric_enzymes]
        if kmod.allosteric_enzymes is not None
        else []
    )
    metabolite_ids = [m.id for m in kmod.metabolites]
    phos_modifying_enzymes = (
        [p.modifying_enzyme_id for p in kmod.phosphorylations]
        if kmod.phosphorylations is not None
        else []
    )
    drain_ids = [
        d.id for d in kmod.reactions if d.mechanism == ReactionMechanism.drain
    ]
    unbalanced_mic_ids, unbalanced_mic_mets, unbalanced_mic_cpts = map(
        list,
        zip(
            *[
                [m.id, m.metabolite_id, m.compartment_id]
                for m in kmod.mics
                if not m.balanced
            ]
        ),
    )
    exp_ids_train = [e.id for e in experiments if e.is_train]
    exp_ids_test = [e.id for e in experiments if e.is_test]
    return ParameterSet(
        dgf=Dgf([metabolite_ids], [[metabolite_ids]], pi.dgf, ii.dgf),
        km=Km([km_ids], [[km_enzs, km_mets, km_cpts]], pi.km, ii.km),
        kcat=Kcat([kcat_ids], [[kcat_enzs, kcat_rxns]], pi.kcat, ii.kcat),
        ki=Ki([ci_ids], [[ci_enzs, ci_rxns, ci_mets, ci_cpts]], pi.ki, ii.ki),
        dissociation_constant=DissociationConstant(
            [dc_ids],
            [[dc_enzs, dc_mets, dc_cpts, dc_mts]],
            pi.dissociation_constant,
            ii.dissociation_constant,
        ),
        transfer_constant=TransferConstant(
            [allosteric_enzyme_ids],
            [[allosteric_enzyme_ids]],
            pi.transfer_constant,
            ii.transfer_constant,
        ),
        kcat_pme=KcatPme(
            [phos_modifying_enzymes],
            [[phos_modifying_enzymes]],
            pi.kcat_pme,
            ii.kcat_pme,
        ),
        drain_train=DrainTrain(
            [exp_ids_train, drain_ids],
            [[exp_ids_train], [drain_ids]],
            pi.drain,
            ii.drain,
        ),
        drain_test=DrainTest(
            [exp_ids_train, drain_ids],
            [[exp_ids_train], [drain_ids]],
            pi.drain,
            ii.drain,
        ),
        conc_enzyme_train=ConcEnzymeTrain(
            [exp_ids_train, enzyme_ids],
            [[exp_ids_train], [enzyme_ids]],
            pi.conc_enzyme,
            ii.conc_enzyme,
            [
                m
                for e in experiments
                for m in e.measurements
                if e.is_train and m.target_type == MeasurementType.ENZYME
            ],
        ),
        conc_enzyme_test=ConcEnzymeTest(
            [exp_ids_train, enzyme_ids],
            [[exp_ids_train], [enzyme_ids]],
            pi.conc_enzyme,
            ii.conc_enzyme,
            [
                m
                for e in experiments
                for m in e.measurements
                if e.is_test and m.target_type == MeasurementType.ENZYME
            ],
        ),
        conc_unbalanced_train=ConcUnbalancedTrain(
            [exp_ids_train, unbalanced_mic_ids],
            [
                [exp_ids_train],
                [unbalanced_mic_mets, unbalanced_mic_cpts],
            ],
            pi.conc_unbalanced,
            ii.conc_unbalanced,
            [
                m
                for e in experiments
                for m in e.measurements
                if e.is_train and m.target_type == MeasurementType.MIC
            ],
        ),
        conc_unbalanced_test=ConcUnbalancedTest(
            [exp_ids_train, unbalanced_mic_ids],
            [
                [exp_ids_train],
                [unbalanced_mic_mets, unbalanced_mic_cpts],
            ],
            pi.conc_unbalanced,
            ii.conc_unbalanced,
            [
                m
                for e in experiments
                for m in e.measurements
                if e.is_test and m.target_type == MeasurementType.MIC
            ],
        ),
        conc_pme_train=ConcPmeTrain(
            [exp_ids_train, phos_modifying_enzymes],
            [[exp_ids_train], [phos_modifying_enzymes]],
            pi.conc_pme,
            ii.conc_pme,
        ),
        conc_pme_test=ConcPmeTest(
            [exp_ids_train, phos_modifying_enzymes],
            [[exp_ids_train], [phos_modifying_enzymes]],
            pi.conc_pme,
            ii.conc_pme,
        ),
        psi_train=PsiTrain([exp_ids_train], [[exp_ids_test]], pi.psi, ii.psi),
        psi_test=PsiTest([exp_ids_train], [[exp_ids_train]], pi.psi, ii.psi),
    )
