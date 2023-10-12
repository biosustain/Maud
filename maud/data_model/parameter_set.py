"""Provides the ParameterSet model.

This is where logic for constructing MaudParameter objects should live.

"""

from pydantic import BaseModel, computed_field

import maud.data_model.parameters as ps
from maud.data_model.experiment import Experiment, MeasurementType
from maud.data_model.hardcoding import ID_SEPARATOR
from maud.data_model.kinetic_model import KineticModel, ReactionMechanism
from maud.data_model.maud_init import InitInput
from maud.data_model.parameter_input import ParameterSetInput


class ParameterSet(BaseModel):
    """the parameters of a maud input."""

    kinetic_model: KineticModel
    experiments: list[Experiment]
    parameter_set_input: ParameterSetInput
    init_input: InitInput

    @computed_field
    def dgf(self) -> ps.Dgf:
        """Add the dgf field."""
        metabolite_ids = [m.id for m in self.kinetic_model.metabolites]
        return ps.Dgf(
            ids=[metabolite_ids],
            split_ids=[[metabolite_ids]],
            user_input=self.parameter_set_input.dgf,
            init_input=self.init_input.dgf,
        )

    @computed_field
    def km(self) -> ps.Km:
        """Add the km field."""
        ids = []
        enzs = []
        mets = []
        cpts = []
        for er in self.kinetic_model.ers:
            rxn = [
                r
                for r in self.kinetic_model.reactions
                if r.id == er.reaction_id
            ][0]
            enz = [
                e for e in self.kinetic_model.enzymes if e.id == er.enzyme_id
            ][0]
            mic_ids = (
                list(rxn.stoichiometry.keys())
                if rxn.mechanism
                != ReactionMechanism.irreversible_michaelis_menten
                else [k for k, v in rxn.stoichiometry.items() if v < 0]
            )
            for mic_id in mic_ids:
                id = ID_SEPARATOR.join([enz.id, mic_id])
                met_id, cpt_id = mic_id.split(ID_SEPARATOR)
                if id not in ids:
                    ids.append(id)
                    enzs.append(enz.id)
                    mets.append(met_id)
                    cpts.append(cpt_id)
        return ps.Km(
            ids=[ids],
            split_ids=[[enzs, mets, cpts]],
            user_input=self.parameter_set_input.km,
            init_input=self.init_input.km,
        )

    @computed_field
    def ki(self) -> ps.Ki:
        """Add the ki field."""
        ids = []
        enzs = []
        rxns = []
        mets = []
        cpts = []
        if self.kinetic_model.competitive_inhibitions is not None:
            for ci in self.kinetic_model.competitive_inhibitions:
                ids.append(ci.id)
                enzs.append(ci.enzyme_id)
                rxns.append(ci.reaction_id)
                mets.append(ci.metabolite_id)
                cpts.append(ci.compartment_id)
        return ps.Ki(
            ids=[ids],
            split_ids=[[enzs, rxns, mets, cpts]],
            user_input=self.parameter_set_input.ki,
            init_input=self.init_input.ki,
        )

    @computed_field
    def kcat(self) -> ps.Kcat:
        """Add the kcat field."""
        ids = []
        enzs = []
        rxns = []
        for er in self.kinetic_model.ers:
            ids.append(er.id)
            enzs.append(er.enzyme_id)
            rxns.append(er.reaction_id)
        return ps.Kcat(
            ids=[ids],
            split_ids=[[enzs, rxns]],
            user_input=self.parameter_set_input.kcat,
            init_input=self.init_input.kcat,
        )

    @computed_field
    def dissociation_constant(self) -> ps.DissociationConstant:
        """Add the dissociation_constant field."""
        ids = []
        enzs = []
        mets = []
        cpts = []
        mts = []
        if self.kinetic_model.allosteries is not None:
            for a in self.kinetic_model.allosteries:
                ids.append(a.id)
                enzs.append(a.enzyme_id)
                mets.append(a.metabolite_id)
                cpts.append(a.compartment_id)
                mts.append(a.modification_type.name)
        return ps.DissociationConstant(
            ids=[ids],
            split_ids=[[enzs, mets, cpts, mts]],
            user_input=self.parameter_set_input.dissociation_constant,
            init_input=self.init_input.dissociation_constant,
        )

    @computed_field
    def transfer_constant(self) -> ps.TransferConstant:
        """Add the transfer_constant field."""
        allosteric_enzyme_ids = (
            [e.id for e in self.kinetic_model.allosteric_enzymes]
            if self.kinetic_model.allosteric_enzymes is not None
            else []
        )
        return ps.TransferConstant(
            ids=[allosteric_enzyme_ids],
            split_ids=[[allosteric_enzyme_ids]],
            user_input=self.parameter_set_input.transfer_constant,
            init_input=self.init_input.transfer_constant,
        )

    @computed_field
    def kcat_pme(self) -> ps.KcatPme:
        """Add the kcat_pme field."""
        phos_modifying_enzymes = (
            [p.modifying_enzyme_id for p in self.kinetic_model.phosphorylations]
            if self.kinetic_model.phosphorylations is not None
            else []
        )
        return ps.KcatPme(
            ids=[phos_modifying_enzymes],
            split_ids=[[phos_modifying_enzymes]],
            user_input=self.parameter_set_input.kcat_pme,
            init_input=self.init_input.kcat_pme,
        )

    @computed_field
    def drain_train(self) -> ps.DrainTrain:
        """Add the drain_train field."""
        drain_ids = [
            d.id
            for d in self.kinetic_model.reactions
            if d.mechanism == ReactionMechanism.drain
        ]
        exp_ids = [e.id for e in self.experiments if e.is_train]
        return ps.DrainTrain(
            ids=[exp_ids, drain_ids],
            split_ids=[[exp_ids], [drain_ids]],
            user_input=self.parameter_set_input.drain,
            init_input=self.init_input.drain,
        )

    @computed_field
    def drain_test(self) -> ps.DrainTest:
        """Add the drain_test field."""
        drain_ids = [
            d.id
            for d in self.kinetic_model.reactions
            if d.mechanism == ReactionMechanism.drain
        ]
        exp_ids = [e.id for e in self.experiments if e.is_test]
        return ps.DrainTest(
            ids=[exp_ids, drain_ids],
            split_ids=[[exp_ids], [drain_ids]],
            user_input=self.parameter_set_input.drain,
            init_input=self.init_input.drain,
        )

    @computed_field
    def conc_enzyme_train(self) -> ps.ConcEnzymeTrain:
        """Add the conc_enzyme_train field."""
        enzyme_ids = [e.id for e in self.kinetic_model.enzymes]
        exp_ids = [e.id for e in self.experiments if e.is_train]
        measurements = [
            m
            for e in self.experiments
            for m in e.measurements
            if e.is_train and m.target_type == MeasurementType.ENZYME
        ]
        return ps.ConcEnzymeTrain(
            ids=[exp_ids, enzyme_ids],
            split_ids=[[exp_ids], [enzyme_ids]],
            user_input=self.parameter_set_input.conc_enzyme,
            init_input=self.init_input.conc_enzyme,
            measurements=measurements,
        )

    @computed_field
    def conc_enzyme_test(self) -> ps.ConcEnzymeTest:
        """Add the conc_enzyme_test field."""
        enzyme_ids = [e.id for e in self.kinetic_model.enzymes]
        exp_ids = [e.id for e in self.experiments if e.is_test]
        measurements = [
            m
            for e in self.experiments
            for m in e.measurements
            if e.is_test and m.target_type == MeasurementType.ENZYME
        ]
        return ps.ConcEnzymeTest(
            ids=[exp_ids, enzyme_ids],
            split_ids=[[exp_ids], [enzyme_ids]],
            user_input=self.parameter_set_input.conc_enzyme,
            init_input=self.init_input.conc_enzyme,
            measurements=measurements,
        )

    @computed_field
    def conc_unbalanced_train(self) -> ps.ConcUnbalancedTrain:
        """Add the conc_unbalanced_train field."""
        exp_ids = [e.id for e in self.experiments if e.is_train]
        measurements = [
            m
            for e in self.experiments
            for m in e.measurements
            if e.is_train and m.target_type == MeasurementType.MIC
        ]
        unbalanced_mic_ids, unbalanced_mic_mets, unbalanced_mic_cpts = map(
            list,
            zip(
                *[
                    [m.id, m.metabolite_id, m.compartment_id]
                    for m in self.kinetic_model.mics
                    if not m.balanced
                ]
            ),
        )
        return ps.ConcUnbalancedTrain(
            ids=[exp_ids, unbalanced_mic_ids],
            split_ids=[
                [exp_ids],
                [unbalanced_mic_mets, unbalanced_mic_cpts],
            ],
            user_input=self.parameter_set_input.conc_unbalanced,
            init_input=self.init_input.conc_unbalanced,
            measurements=measurements,
        )

    @computed_field
    def conc_unbalanced_test(self) -> ps.ConcUnbalancedTest:
        """Add the conc_unbalanced_test field."""
        exp_ids = [e.id for e in self.experiments if e.is_test]
        measurements = [
            m
            for e in self.experiments
            for m in e.measurements
            if e.is_test and m.target_type == MeasurementType.MIC
        ]
        unbalanced_mic_ids, unbalanced_mic_mets, unbalanced_mic_cpts = map(
            list,
            zip(
                *[
                    [m.id, m.metabolite_id, m.compartment_id]
                    for m in self.kinetic_model.mics
                    if not m.balanced
                ]
            ),
        )
        return ps.ConcUnbalancedTest(
            ids=[exp_ids, unbalanced_mic_ids],
            split_ids=[
                [exp_ids],
                [unbalanced_mic_mets, unbalanced_mic_cpts],
            ],
            user_input=self.parameter_set_input.conc_unbalanced,
            init_input=self.init_input.conc_unbalanced,
            measurements=measurements,
        )

    @computed_field
    def conc_pme_train(self) -> ps.ConcPmeTrain:
        """Add the conc_pme_train field."""
        exp_ids = [e.id for e in self.experiments if e.is_train]
        pme_ids = (
            [p.modifying_enzyme_id for p in self.kinetic_model.phosphorylations]
            if self.kinetic_model.phosphorylations is not None
            else []
        )
        return ps.ConcPmeTrain(
            ids=[exp_ids, pme_ids],
            split_ids=[[exp_ids], [pme_ids]],
            user_input=self.parameter_set_input.conc_pme,
            init_input=self.init_input.conc_pme,
        )

    @computed_field
    def conc_pme_test(self) -> ps.ConcPmeTest:
        """Add the conc_pme_test field."""
        exp_ids = [e.id for e in self.experiments if e.is_test]
        pme_ids = (
            [p.modifying_enzyme_id for p in self.kinetic_model.phosphorylations]
            if self.kinetic_model.phosphorylations is not None
            else []
        )
        return ps.ConcPmeTest(
            ids=[exp_ids, pme_ids],
            split_ids=[[exp_ids], [pme_ids]],
            user_input=self.parameter_set_input.conc_pme,
            init_input=self.init_input.conc_pme,
        )

    @computed_field
    def psi_train(self) -> ps.PsiTrain:
        """Add the psi_train field."""
        exp_ids = [e.id for e in self.experiments if e.is_train]
        return ps.PsiTrain(
            ids=[exp_ids],
            split_ids=[[exp_ids]],
            user_input=self.parameter_set_input.psi,
            init_input=self.init_input.psi,
        )

    @computed_field
    def psi_test(self) -> ps.PsiTest:
        """Add the psi_test field."""
        exp_ids = [e.id for e in self.experiments if e.is_test]
        return ps.PsiTest(
            ids=[exp_ids],
            split_ids=[[exp_ids]],
            user_input=self.parameter_set_input.psi,
            init_input=self.init_input.psi,
        )
