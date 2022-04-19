from dataclasses import field, fields
from math import isnan
from typing import Dict, Sequence, Union

from pydantic import StrictFloat, StrictInt, root_validator
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
    name: str
    value: StanDataValue

    @root_validator
    def no_nans_allowed(cls, values):
        msg = f"StanData {values['name']} has non-numbers: {values['value']}"
        if isinstance(values["value"], Sequence):
            flat = recursively_flatten_list(values["value"])
            assert not any([isnan(x) for x in flat]), msg
        else:
            assert not isnan(values["value"]), msg
        return values


@dataclass
class StanInputTrain:
    # priors
    prior_loc_dgf: StanData
    prior_cov_dgf: StanData
    priors_kcat: StanData
    priors_km: StanData
    priors_ki: StanData
    priors_dissociation_constant: StanData
    priors_transfer_constant: StanData
    priors_kcat_phos: StanData
    priors_conc_phos: StanData
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
    edge_type: StanData
    edge_to_enzyme: StanData
    edge_to_tc: StanData
    edge_to_drain: StanData
    edge_to_reaction: StanData
    water_stoichiometry: StanData
    mic_to_met: StanData
    subunits: StanData
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
    enzyme_knockout_long: StanData
    enzyme_knockout_bounds: StanData
    phosphorylation_knockout_long: StanData
    phosphorylation_knockout_bounds: StanData

    # user configuration
    conc_init: StanData
    likelihood: StanData
    reject_non_steady: StanData
    steady_state_threshold_abs: StanData
    steady_state_threshold_rel: StanData
    rel_tol: StanData
    abs_tol: StanData
    timepoint: StanData
    max_num_steps: StanData

    # Data that can be inferred
    N_mic: StanData = field(init=False)
    N_edge_sub: StanData = field(init=False)
    N_edge_prod: StanData = field(init=False)
    N_edge: StanData = field(init=False)
    N_unbalanced: StanData = field(init=False)
    N_metabolite: StanData = field(init=False)
    N_km: StanData = field(init=False)
    N_enzyme: StanData = field(init=False)
    N_phosphorylation: StanData = field(init=False)
    N_experiment: StanData = field(init=False)
    N_competitive_inhibition: StanData = field(init=False)
    N_allostery: StanData = field(init=False)
    N_allosteric_enzyme: StanData = field(init=False)
    N_drain: StanData = field(init=False)
    N_flux_measurement: StanData = field(init=False)
    N_enzyme_measurement: StanData = field(init=False)
    N_conc_measurement: StanData = field(init=False)
    N_enzyme_knockout: StanData = field(init=False)
    N_phosphorylation_knockout: StanData = field(init=False)

    # final output
    stan_input_dict: StanInputDict = field(init=False)

    def __post_init__(self):
        self.N_mic = StanData(name="N_mic", value=len(self.S.value))
        self.N_edge_sub = StanData(
            name="N_edge_sub", value=len(self.sub_by_edge_long.value)
        )
        self.N_edge_prod = StanData(
            name="N_edge_sub", value=len(self.prod_by_edge_long.value)
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
        self.N_enzyme = StanData(
            name="N_enzyme", value=len(self.priors_kcat.value[0])
        )
        self.N_phosphorylation = StanData(
            name="N_phosphorylation", value=len(self.priors_kcat_phos.value[0])
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
        self.N_phosphorylation_knockout = StanData(
            name="N_phosphorylation_knockout",
            value=len(self.phosphorylation_knockout_long.value),
        )
        self.stan_input_dict = {
            f.name: getattr(self, f.name).value
            for f in fields(self)
            if f.name != "stan_input_dict"
        }


@dataclass
class StanInputTest:
    # priors
    priors_conc_phos: StanData
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
    edge_type: StanData
    edge_to_enzyme: StanData
    edge_to_tc: StanData
    edge_to_drain: StanData
    edge_to_reaction: StanData
    water_stoichiometry: StanData
    mic_to_met: StanData
    subunits: StanData
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
    enzyme_knockout_long: StanData
    enzyme_knockout_bounds: StanData
    phosphorylation_knockout_long: StanData
    phosphorylation_knockout_bounds: StanData

    # user configuration
    conc_init: StanData
    likelihood: StanData
    rel_tol: StanData
    abs_tol: StanData
    timepoint: StanData
    max_num_steps: StanData

    # Data that can be inferred
    N_mic: StanData = field(init=False)
    N_edge_sub: StanData = field(init=False)
    N_edge_prod: StanData = field(init=False)
    N_edge: StanData = field(init=False)
    N_unbalanced: StanData = field(init=False)
    N_experiment: StanData = field(init=False)
    N_competitive_inhibition: StanData = field(init=False)
    N_allostery: StanData = field(init=False)
    N_phosphorylation: StanData = field(init=False)
    N_drain: StanData = field(init=False)
    N_enzyme: StanData = field(init=False)
    N_enzyme_knockout: StanData = field(init=False)
    N_phosphorylation_knockout: StanData = field(init=False)

    # final output
    stan_input_dict: StanInputDict = field(init=False)

    def __post_init__(self):
        self.N_mic = StanData(name="N_mic", value=len(self.S.value))
        self.N_edge_sub = StanData(
            name="N_edge_sub", value=len(self.sub_by_edge_long.value)
        )
        self.N_edge_prod = StanData(
            name="N_edge_sub", value=len(self.prod_by_edge_long.value)
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
        self.N_phosphorylation_knockout = StanData(
            name="N_phosphorylation_knockout",
            value=len(self.phosphorylation_knockout_long.value),
        )
        self.stan_input_dict = {
            f.name: getattr(self, f.name).value
            for f in fields(self)
            if f.name != "stan_input_dict"
        }
