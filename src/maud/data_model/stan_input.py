from dataclasses import field
from typing import Union, Sequence
from pydantic.dataclasses import dataclass


StanDataValueBase = Union[int, float]
StanDataValue = Union[
    StanDataValueBase,
    Sequence[StanDataValueBase],
    Sequence[Sequence[StanDataValueBase]],
    Sequence[Sequence[Sequence[StanDataValueBase]]],
]


@dataclass
class StanData:
    value: StanDataValue


@dataclass
class StanInput:
    # priors
    prior_loc_dgf: StanData
    prior_cov_dgf: StanData
    priors_kcat: StanData
    priors_km: StanData
    priors_ki: StanData
    priors_dc: StanData
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

    # Data that can be inferred
    N_mic: StanData = field(init=False)
    N_edge: StanData = field(init=False)
    N_unbalanced: StanData = field(init=False)
    N_metabolite: StanData = field(init=False)
    N_km: StanData = field(init=False)
    N_reaction: StanData = field(init=False)
    N_enzyme: StanData = field(init=False)
    N_phosphorylation_enzymes: StanData = field(init=False)
    N_experiment: StanData = field(init=False)
    N_ci: StanData = field(init=False)
    N_ai: StanData = field(init=False)
    N_aa: StanData = field(init=False)
    N_ae: StanData = field(init=False)
    N_drain: StanData = field(init=False)
    N_flux_measurement: StanData = field(init=False)
    N_enzyme_measurement: StanData = field(init=False)
    N_conc_measurement: StanData = field(init=False)

    def __post_init__(self):
        self.N_mic = StanData(len(self.S.value))
        self.N_edge = StanData(len(self.S.value[0]))
        self.N_unbalanced = StanData(len(self.priors_conc_unbalanced.value[0][0]))
        self.N_metabolite = StanData(len(self.prior_loc_dgf.value))
        self.N_km = StanData(len(self.priors_km.value[0]))
        self.N_reaction = StanData(len(self.priors_kcat.value[0] + len(self.priors_drain.value[0])))
        self.N_enzyme = StanData(len(self.priors_kcat.value[0]))
        self.N_phosphorylation_enzymes = StanData(len(self.priors_kcat_phos.value[0]))
        self.N_experiment = StanData(len(self.priors_conc_unbalanced.value[0]))
        self.N_ci = StanData(len(self.priors_ki.value[0]))
        self.N_ai = StanData(name="N_ai", value=len(self.diss_t.ids.value[0]))
        self.N_aa = StanData(name="N_aa", value=len(self.diss_r.ids.value[0]))
        self.N_ae = StanData(len(self.priors_transfer_constant.value[0]))
        self.N_drain = StanData(len(self.priors_drain.value[0][0]))
        self.N_flux_measurement = StanData(len(self.yflux.value[0]))
        self.N_enzyme_measurement = StanData(len(self.yenz.value[0]))
        self.N_conc_measurement = StanData(len(self.yconc.value[0]))
