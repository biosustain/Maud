"""Definitions of inputs for Maud's Stan models."""

from dataclasses import fields
from math import isnan
from typing import Dict, Sequence, Union

from pydantic import Field, StrictFloat, StrictInt, root_validator
from pydantic.dataclasses import dataclass

from maud.utils import recursively_flatten_list


# Strict types are used here to stop pydantic from doing unwanted conversions.
# see these links for more discussion:
# - https://pydantic-docs.helpmanual.io/usage/models/#data-conversion
# - https://pydantic-docs.helpmanual.io/usage/types/#strict-types
StanDataValueBase = Union[StrictFloat, StrictInt]
StanDataValue = Union[
    StanDataValueBase,
    Sequence[StanDataValueBase],
    Sequence[Sequence[StanDataValueBase]],
    Sequence[Sequence[Sequence[StanDataValueBase]]],
]
StanInputDict = Dict[str, StanDataValue]


@dataclass
class StanData:
    """Something that will go in Stan's data block, in Python form."""

    name: str
    value: StanDataValue

    @root_validator
    def no_nans_allowed(cls, values):
        """Check that Stan input data isn't null and doesn't contain nulls."""
        msg = f"StanData {values['name']} has non-numbers: {values['value']}"
        if isinstance(values["value"], Sequence):
            flat = recursively_flatten_list(values["value"])
            assert not any([isnan(x) for x in flat]), msg
        else:
            assert not isnan(values["value"]), msg
        return values


@dataclass
class StanInputTrain:
    """Input for model.stan, i.e. the training model."""

    # priors
    prior_loc_dgf: StanData
    prior_cov_dgf: StanData
    priors_kcat: StanData
    priors_km: StanData
    priors_ki: StanData
    priors_dissociation_constant: StanData
    priors_transfer_constant: StanData
    priors_kcat_pme: StanData
    priors_psi: StanData
    priors_conc_pme: StanData
    priors_conc_unbalanced: StanData
    priors_conc_enzyme: StanData
    priors_drain: StanData

    # measurements
    yconc: StanData
    sigma_conc: StanData
    experiment_yconc: StanData
    mic_ix_yconc: StanData
    yflux: StanData
    sigma_flux: StanData
    experiment_yflux: StanData
    reaction_yflux: StanData
    yenz: StanData
    sigma_enz: StanData
    experiment_yenz: StanData
    enzyme_yenz: StanData

    # network properties
    S: StanData
    N_reaction: StanData
    balanced_mic_ix: StanData
    unbalanced_mic_ix: StanData
    ci_mic_ix: StanData
    edge_type: StanData
    edge_to_enzyme: StanData
    edge_to_tc: StanData
    edge_to_drain: StanData
    edge_to_reaction: StanData
    water_stoichiometry: StanData
    transported_charge: StanData
    mic_to_met: StanData
    subunits: StanData
    temperature: StanData
    sub_by_edge_long: StanData
    sub_by_edge_bounds: StanData
    prod_by_edge_long: StanData
    prod_by_edge_bounds: StanData
    sub_km_ix_by_edge_long: StanData
    sub_km_ix_by_edge_bounds: StanData
    prod_km_ix_by_edge_long: StanData
    prod_km_ix_by_edge_bounds: StanData
    ci_ix_long: StanData
    ci_ix_bounds: StanData
    allostery_ix_long: StanData
    allostery_ix_bounds: StanData
    allostery_type: StanData
    allostery_mic: StanData
    phosphorylation_ix_long: StanData
    phosphorylation_ix_bounds: StanData
    phosphorylation_type: StanData
    phosphorylation_pme: StanData
    enzyme_knockout_long: StanData
    enzyme_knockout_bounds: StanData
    pme_knockout_long: StanData
    pme_knockout_bounds: StanData

    # user configuration
    conc_init: StanData
    likelihood: StanData
    drain_small_conc_corrector: StanData
    reject_non_steady: StanData
    steady_state_threshold_abs: StanData
    steady_state_threshold_rel: StanData
    rel_tol: StanData
    abs_tol: StanData
    timepoint: StanData
    max_num_steps: StanData

    # Data that can be inferred
    N_mic: StanData = Field(init=False, exclude=True)
    N_edge_sub: StanData = Field(init=False, exclude=True)
    N_edge_prod: StanData = Field(init=False, exclude=True)
    N_edge: StanData = Field(init=False, exclude=True)
    N_unbalanced: StanData = Field(init=False, exclude=True)
    N_metabolite: StanData = Field(init=False, exclude=True)
    N_km: StanData = Field(init=False, exclude=True)
    N_sub_km: StanData = Field(init=False, exclude=True)
    N_prod_km: StanData = Field(init=False, exclude=True)
    N_enzyme: StanData = Field(init=False, exclude=True)
    N_phosphorylation: StanData = Field(init=False, exclude=True)
    N_pme: StanData = Field(init=False, exclude=True)
    N_experiment: StanData = Field(init=False, exclude=True)
    N_competitive_inhibition: StanData = Field(init=False, exclude=True)
    N_allostery: StanData = Field(init=False, exclude=True)
    N_allosteric_enzyme: StanData = Field(init=False, exclude=True)
    N_drain: StanData = Field(init=False, exclude=True)
    N_flux_measurement: StanData = Field(init=False, exclude=True)
    N_enzyme_measurement: StanData = Field(init=False, exclude=True)
    N_conc_measurement: StanData = Field(init=False, exclude=True)
    N_enzyme_knockout: StanData = Field(init=False, exclude=True)
    N_pme_knockout: StanData = Field(init=False, exclude=True)

    # final output
    stan_input_dict: StanInputDict = Field(init=False, exclude=True)

    def __post_init__(self):
        """Add fields that can be inferred from other ones."""
        self.N_mic = StanData(name="N_mic", value=len(self.S.value))
        self.N_edge_sub = StanData(
            name="N_edge_sub", value=len(self.sub_by_edge_long.value)
        )
        self.N_edge_prod = StanData(
            name="N_edge_prod", value=len(self.prod_by_edge_long.value)
        )
        self.N_edge = StanData(name="N_edge", value=len(self.S.value[0]))
        self.N_unbalanced = StanData(
            name="N_unbalanced",
            value=len(self.priors_conc_unbalanced.value[0][0]),
        )
        self.N_metabolite = StanData(
            name="N_metabolite", value=len(self.prior_loc_dgf.value)
        )
        self.N_km = StanData(name="N_km", value=len(self.priors_km.value[0]))
        self.N_sub_km = StanData(
            name="N_sub_km", value=len(self.sub_km_ix_by_edge_long.value)
        )
        self.N_prod_km = StanData(
            name="N_prod_km", value=len(self.prod_km_ix_by_edge_long.value)
        )
        self.N_enzyme = StanData(
            name="N_enzyme", value=len(self.priors_kcat.value[0])
        )
        self.N_phosphorylation = StanData(
            name="N_phosphorylation",
            value=len(self.phosphorylation_type.value),
        )
        self.N_pme = StanData(
            name="N_pme", value=len(self.priors_kcat_pme.value[0])
        )
        self.N_experiment = StanData(
            name="N_experiment",
            value=len(self.priors_conc_unbalanced.value[0]),
        )
        self.N_competitive_inhibition = StanData(
            name="N_competitive_inhibition", value=len(self.priors_ki.value[0])
        )
        self.N_allostery = StanData(
            name="N_allostery", value=len(self.allostery_type.value)
        )
        self.N_allosteric_enzyme = StanData(
            name="N_allosteric_enzyme",
            value=len(self.priors_transfer_constant.value[0]),
        )
        self.N_drain = StanData(
            name="N_drain", value=len(self.priors_drain.value[0][0])
        )
        self.N_flux_measurement = StanData(
            name="N_flux_measurement", value=len(self.yflux.value)
        )
        self.N_enzyme_measurement = StanData(
            name="N_enzyme_measurement", value=len(self.yenz.value)
        )
        self.N_conc_measurement = StanData(
            name="N_conc_measurement", value=len(self.yconc.value)
        )
        self.N_enzyme_knockout = StanData(
            name="N_enzyme_knockout",
            value=len(self.enzyme_knockout_long.value),
        )
        self.N_pme_knockout = StanData(
            name="N_pme_knockout",
            value=len(self.pme_knockout_long.value),
        )
        self.stan_input_dict = {
            f.name: getattr(self, f.name).value
            for f in fields(self)
            if f.name != "stan_input_dict"
        }


@dataclass
class StanInputTest:
    """Input for out_of_sample_model.stan, i.e. the test model."""

    # priors
    priors_conc_pme: StanData
    priors_conc_unbalanced: StanData
    priors_conc_enzyme: StanData
    priors_drain: StanData

    # network properties
    S: StanData
    N_metabolite: StanData
    N_km: StanData
    N_reaction: StanData
    N_allosteric_enzyme: StanData
    balanced_mic_ix: StanData
    unbalanced_mic_ix: StanData
    ci_mic_ix: StanData
    edge_type: StanData
    edge_to_enzyme: StanData
    edge_to_tc: StanData
    edge_to_drain: StanData
    edge_to_reaction: StanData
    water_stoichiometry: StanData
    mic_to_met: StanData
    subunits: StanData
    temperature: StanData
    sub_by_edge_long: StanData
    sub_by_edge_bounds: StanData
    prod_by_edge_long: StanData
    prod_by_edge_bounds: StanData
    sub_km_ix_by_edge_long: StanData
    sub_km_ix_by_edge_bounds: StanData
    prod_km_ix_by_edge_long: StanData
    prod_km_ix_by_edge_bounds: StanData
    ci_ix_long: StanData
    ci_ix_bounds: StanData
    allostery_ix_long: StanData
    allostery_ix_bounds: StanData
    allostery_type: StanData
    allostery_mic: StanData
    phosphorylation_ix_long: StanData
    phosphorylation_ix_bounds: StanData
    phosphorylation_type: StanData
    phosphorylation_pme: StanData
    enzyme_knockout_long: StanData
    enzyme_knockout_bounds: StanData
    pme_knockout_long: StanData
    pme_knockout_bounds: StanData

    # user configuration
    conc_init: StanData
    likelihood: StanData
    drain_small_conc_corrector: StanData
    rel_tol: StanData
    abs_tol: StanData
    timepoint: StanData
    max_num_steps: StanData

    # Data that can be inferred
    N_mic: StanData = Field(init=False, exclude=True)
    N_edge_sub: StanData = Field(init=False, exclude=True)
    N_edge_prod: StanData = Field(init=False, exclude=True)
    N_edge: StanData = Field(init=False, exclude=True)
    N_unbalanced: StanData = Field(init=False, exclude=True)
    N_experiment: StanData = Field(init=False, exclude=True)
    N_competitive_inhibition: StanData = Field(init=False, exclude=True)
    N_allostery: StanData = Field(init=False, exclude=True)
    N_phosphorylation: StanData = Field(init=False, exclude=True)
    N_pme: StanData = Field(init=False, exclude=True)
    N_drain: StanData = Field(init=False, exclude=True)
    N_enzyme: StanData = Field(init=False, exclude=True)
    N_enzyme_knockout: StanData = Field(init=False, exclude=True)
    N_pme_knockout: StanData = Field(init=False, exclude=True)
    N_sub_km: StanData = Field(init=False, exclude=True)
    N_prod_km: StanData = Field(init=False, exclude=True)

    # final output
    stan_input_dict: StanInputDict = Field(init=False, exclude=True)

    def __post_init__(self):
        """Add fields that can be inferred from other ones."""
        self.N_mic = StanData(name="N_mic", value=len(self.S.value))
        self.N_edge_sub = StanData(
            name="N_edge_sub", value=len(self.sub_by_edge_long.value)
        )
        self.N_edge_prod = StanData(
            name="N_edge_prod", value=len(self.prod_by_edge_long.value)
        )
        self.N_edge = StanData(name="N_edge", value=len(self.S.value[0]))
        self.N_unbalanced = StanData(
            name="N_unbalanced",
            value=len(self.priors_conc_unbalanced.value[0][0]),
        )
        self.N_enzyme = StanData("N_enzyme", len(self.subunits.value))
        self.N_phosphorylation = StanData(
            "N_phosphorylation", len(self.phosphorylation_ix_long.value)
        )
        self.N_pme = StanData(
            name="N_pme", value=len(self.priors_conc_pme.value[0][0])
        )
        self.N_experiment = StanData(
            "N_experiment",
            len(self.priors_conc_unbalanced.value[0]),
        )
        self.N_competitive_inhibition = StanData(
            name="N_competitive_inhibition", value=len(self.ci_ix_long.value)
        )
        self.N_allostery = StanData(
            name="N_allostery", value=len(self.allostery_ix_long.value)
        )
        self.N_drain = StanData(
            name="N_drain", value=len(self.priors_drain.value[0][0])
        )
        self.N_enzyme_knockout = StanData(
            name="N_enzyme_knockout",
            value=len(self.enzyme_knockout_long.value),
        )
        self.N_pme_knockout = StanData(
            name="N_pme_knockout", value=len(self.pme_knockout_long.value)
        )
        self.N_sub_km = StanData(
            name="N_sub_km", value=len(self.sub_km_ix_by_edge_long.value)
        )
        self.N_prod_km = StanData(
            name="N_prod_km", value=len(self.prod_km_ix_by_edge_long.value)
        )
        self.stan_input_dict = {
            f.name: getattr(self, f.name).value
            for f in fields(self)
            if f.name != "stan_input_dict"
        }
