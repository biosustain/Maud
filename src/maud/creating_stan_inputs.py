from typing import Dict, Iterable, List, Tuple, Union

import numpy as np
import pandas as pd

from maud.data_model.hardcoding import ID_SEPARATOR
from maud.data_model.kinetic_model import (
    Allostery,
    CompetitiveInhibition,
    Enzyme,
    EnzymeReaction,
    KineticModel,
    Metabolite,
    MetaboliteInCompartment,
    Phosphorylation,
    Reaction,
)
from maud.data_model.maud_config import MaudConfig
from maud.data_model.measurement_set import Experiment, MeasurementSet
from maud.data_model.prior_set import IndPrior1d, IndPrior2d, PriorSet
from maud.data_model.stan_input import StanData, StanInput
from maud.data_model.stan_variable_set import DissociationConstant, Ki, Km


CodifiableMaudObject = Union[
    Allostery,
    MetaboliteInCompartment,
    Metabolite,
    EnzymeReaction,
    Enzyme,
    Phosphorylation,
    Reaction,
    Experiment,
    CompetitiveInhibition
]


def codify_maud_object(lx: Iterable[CodifiableMaudObject]) -> Dict[str, int]:
    """Turn an iterable of Maud objects into a dictionary of Stan codes."""
    return dict(zip([x.id for x in lx], range(1, len(list(lx)) + 1)))


def create_stan_input(
    ms: MeasurementSet,
    priors: PriorSet,
    kinetic_model: KineticModel,
    config: MaudConfig
) -> StanInput:
    """Get the input to inference_model.stan from a MaudInput object.

    :param mi: a MaudInput object

    """
    S = kinetic_model.stoichiometric_matrix
    mics = kinetic_model.mics
    edges = kinetic_model.edges
    allosteries = kinetic_model.allosteries
    enzymes = kinetic_model.enzymes
    reactions = kinetic_model.reactions
    reaction_by_id = {r.id: r for r in reactions}
    metabolite_codes = codify_maud_object(kinetic_model.metabolites)
    mic_codes = codify_maud_object(mics)
    enzyme_codes = codify_maud_object(enzymes)
    experiment_codes = codify_maud_object(ms.experiments)
    # phos_enz_codes
    if allosteries is not None:
        tc_codes = codify_maud_object([
            e for e in enzymes if any(a.enzyme_id == e.id for a in allosteries)
        ])
        allostery_codes = codify_maud_object(allosteries)
    else:
        tc_codes = dict()
        allostery_codes = dict()
    if kinetic_model.competitive_inhibitions is not None:
        ci_codes = codify_maud_object(kinetic_model.competitive_inhibitions)
    else:
        ci_codes = dict()
    if kinetic_model.phosphorylations is not None:
        phosphorylation_codes = codify_maud_object(
            kinetic_model.phosphorylations)
    else:
        phosphorylation_codes = dict()
    drain_codes = codify_maud_object(kinetic_model.drains)
    reaction_codes = codify_maud_object(reactions)
    reaction_to_mechanism = {rxn.id: rxn.mechanism for rxn in reactions}
    edge_mechanism_code = [
        reaction_to_mechanism[e.id].value if isinstance(e, Reaction)
        else reaction_to_mechanism[e.reaction_id].value
        for e in edges
    ]
    edge_enzyme_code = [
        enzyme_codes[e.id] if isinstance(e, Enzyme) else 0
        for e in edges
    ]
    edge_drain_code = [
        drain_codes[e.id] if isinstance(e, Reaction) else 0
        for e in edges
    ]
    edge_reaction_code = [
        reaction_codes[e.id] if isinstance(e, Reaction)
        else reaction_codes[e.reaction_id]
        for e in edges
    ]
    edge_water_stoichiometry = [
        reaction_by_id[e.reaction_id].water_stoichiometry
        if isinstance(e, EnzymeReaction) else 0 for e in edges
    ]
    mic_met_code = [metabolite_codes[m.metabolite_id] for m in mics]
    balanced_mic_codes = [mic_codes[mic.id] for mic in mics if mic.balanced]
    unbalanced_mic_codes = [mic_codes[mic.id]
                            for mic in mics if not mic.balanced]
    edge_tc_code = [
        tc_codes[e.enzyme_id]
        if isinstance(e, EnzymeReaction) and e.enzyme_id in tc_codes.keys()
        else 0 for e in edges
    ]
    # ragged arrays
    sub_by_edge = [list(row.loc[row < 0].index) for _, row in S.T.iterrows()]
    prod_by_edge = [list(row.loc[row < 0].index) for _, row in S.T.iterrows()]
    enz_ko_by_experiment = [
        [
            enzyme_codes[eko.enzyme_id]
            for eko in ms.enzyme_knockouts
            if eko.experiment_id == e
        ] for e in ms.experiments
    ] if ms.enzyme_knockouts is not None else [[] for e in ms.experiments]
    phos_ko_by_experiment = [
        [
            enzyme_codes[pko.enzyme_id]
            for pko in ms.phosphorylation_knockouts
            if pko.experiment_id == e
        ] for e in ms.experiments
    ] if ms.phosphorylation_knockouts is not None else [
        [] for e in ms.experiments
    ]
    ci_by_edge: List[List[int]] = [[] for er in kinetic_model.edges]
    if kinetic_model.competitive_inhibitions is not None:
        for ci in kinetic_model.competitive_inhibitions:
            edge_pos = next(
                i for i, e in enumerate(kinetic_model.edges)
                if (
                    isinstance(e, EnzymeReaction)
                    and e.enzyme_id == ci.enzyme_id
                )
            )
            ci_by_edge[edge_pos].append(ci_codes[ci.id])
    allostery_by_edge: List[List[int]] = [[] for er in kinetic_model.edges]
    if kinetic_model.allosteries is not None:
        allostery_type = [
            a.modification_type.value for a in kinetic_model.allosteries]
        allostery_mic = [mic_codes[a.mic_id]
                         for a in kinetic_model.allosteries]
        for a in kinetic_model.allosteries:
            edge_pos = next(  # is this right???
                i for i, er in enumerate(kinetic_model.edges)
                if (
                    isinstance(er, EnzymeReaction)
                    and er.enzyme_id == a.enzyme_id
                )
            )
            allostery_by_edge[edge_pos].append(allostery_codes[a.id])
    else:
        allostery_type = []
        allostery_mic = []
    phosphorylation_by_edge: List[List[int]] = [[]
                                                for er in kinetic_model.edges]
    if kinetic_model.phosphorylations is not None:
        phosphorylation_type = [
            p.modification_type.value for p in kinetic_model.phosphorylations]
        for p in kinetic_model.phosphorylations:
            edge_pos = next(
                i for i, er in enumerate(kinetic_model.edges)
                if (
                    isinstance(er, EnzymeReaction)
                    and er.enzyme_id == p.enzyme_id
                )
            )
            phosphorylation_by_edge[edge_pos].append(
                phosphorylation_codes[p.id])
    else:
        phosphorylation_type = []
    sub_by_edge_long, sub_by_edge_bounds = encode_ragged(sub_by_edge)
    prod_by_edge_long, prod_by_edge_bounds = encode_ragged(prod_by_edge)
    enz_knockout_long, enz_knockout_bounds = encode_ragged(
        enz_ko_by_experiment)
    phos_knockout_long, phos_knockout_bounds = encode_ragged(
        phos_ko_by_experiment)
    ci_ix_long, ci_ix_bounds = encode_ragged(ci_by_edge)
    allostery_ix_long, allostery_ix_bounds = encode_ragged(allostery_by_edge)
    phosphorylation_ix_long, phosphorylation_ix_bounds = encode_ragged(
        phosphorylation_by_edge)
    return StanInput(
        # priors
        prior_loc_dgf=StanData(priors.dgf.location.tolist()),
        prior_cov_dgf=StanData(priors.dgf.covariance_matrix.tolist()),
        priors_kcat=StanData(unpack_priors(priors.kcat)),
        priors_km=StanData(unpack_priors(priors.km)),
        priors_ki=StanData(unpack_priors(priors.ki)),
        priors_dc=StanData(unpack_priors(priors.dissociation_constant)),
        priors_kcat_phos=StanData(unpack_priors(priors.kcat_phos)),
        priors_transfer_constant=StanData(
            unpack_priors(priors.transfer_constant)),
        priors_conc_unbalanced=StanData(unpack_priors(priors.conc_unbalanced)),
        priors_conc_enzyme=StanData(unpack_priors(priors.conc_enzyme)),
        priors_conc_phos=StanData(unpack_priors(priors.conc_phos)),
        priors_drain=StanData(unpack_priors(priors.drain)),
        # measurements
        yconc=StanData(ms.yconc["measurement"].tolist()),
        sigma_conc=StanData(ms.yconc["error_scale"].tolist()),
        experiment_yconc=StanData(
            ms.yconc["experiment_id"].map(experiment_codes).tolist()),
        mic_ix_yconc=StanData(ms.yconc.reset_index()[
                              "target_id"].map(mic_codes).tolist()),
        yflux=StanData(ms.yflux["measurement"].tolist()),
        sigma_flux=StanData(ms.yflux["error_scale"].tolist()),
        experiment_yflux=StanData(ms.yflux.reset_index(
        )["experiment_id"].map(experiment_codes).tolist()),
        reaction_yflux=StanData(ms.yflux.reset_index(
        )["target_id"].map(reaction_codes).tolist()),
        yenz=StanData(ms.yenz["measurement"].tolist()),
        sigma_enz=StanData(ms.yenz["error_scale"].tolist()),
        experiment_yenz=StanData(ms.yenz.reset_index(
        )["experiment_id"].map(experiment_codes).tolist()),
        enzyme_yenz=StanData(ms.yenz.reset_index()[
                             "target_id"].map(enzyme_codes).tolist()),
        # network properties
        S=StanData(S.tolist()),
        balanced_mic_ix=StanData(balanced_mic_codes),
        unbalanced_mic_ix=StanData(unbalanced_mic_codes),
        edge_type=StanData(edge_mechanism_code),
        edge_to_enzyme=StanData(edge_enzyme_code),
        edge_to_tc=StanData(edge_tc_code),
        edge_to_drain=StanData(edge_drain_code),
        edge_to_reaction=StanData(edge_reaction_code),
        water_stoichiometry=StanData(edge_water_stoichiometry),
        mic_to_met=StanData(mic_met_code),
        subunits=StanData([e.subunits for e in kinetic_model.enzymes]),
        sub_by_edge_long=StanData(sub_by_edge_long),
        sub_by_edge_bounds=StanData(sub_by_edge_bounds),
        prod_by_edge_long=StanData(prod_by_edge_long),
        prod_by_edge_bounds=StanData(prod_by_edge_bounds),
        ci_ix_long=StanData(ci_ix_long),
        ci_ix_bounds=StanData(ci_ix_bounds),
        allostery_ix_long=StanData(allostery_ix_long),
        allostery_ix_bounds=StanData(allostery_ix_bounds),
        allostery_type=StanData(allostery_type),
        allostery_mic=StanData(allostery_mic),
        phosphorylation_ix_long=StanData(phosphorylation_ix_long),
        phosphorylation_ix_bounds=StanData(phosphorylation_ix_bounds),
        phosphorylation_type=StanData(phosphorylation_type),
        enzyme_knockout_long=StanData(enz_knockout_long),
        enzyme_knockout_bounds=StanData(enz_knockout_bounds),
        phosphorylation_knockout_long=StanData(phos_knockout_long),
        phosphorylation_knockout_bounds=StanData(phos_knockout_bounds),
        # configuration
        conc_init=get_conc_init(ms, kinetic_model, priors, config),
        likelihood=StanData(config.likelihood),
        reject_non_steady=StanData(config.reject_non_steady),
        steady_state_threshold_abs=StanData(config.steady_state_threshold_abs),
        steady_state_threshold_rel=StanData(config.steady_state_threshold_rel),
        rel_tol=StanData(config.ode_config.rel_tol),
        abs_tol=StanData(config.ode_config.abs_tol),
        timepoint=StanData(config.ode_config.timepoint),
    )


def encode_ragged(ragged: List[List]) -> Tuple[List, List]:
    """Encode a ragged list of lists in a Stan friendly format.

    Specifically, the return value is a flat list with all the data, and a list
    of bounds with the start and end points (1-indexed) of each entry.

    """
    flat: list = []
    bounds: list = []
    ticker = 1
    for sublist in ragged:
        flat = flat + sublist
        bounds.append([ticker, ticker + len(sublist) - 1])
        ticker += len(sublist)
    return flat, bounds


def get_lookup(k: Union[Km, Ki], S: pd.DataFrame) -> List[int]:
    out = pd.DataFrame(np.zeros(S.shape), index=S.index, columns=S.columns)
    for i, (mic_id, er_id) in enumerate(zip(k.mic_ids, k.er_ids)):
        out.loc[mic_id, er_id] = i + 1
    return out.values.astype(int).tolist()


def get_diss_lookup(d: DissociationConstant, S: pd.DataFrame) -> List[int]:
    out = pd.DataFrame(np.zeros(S.shape), index=S.index, columns=S.columns)
    for i, (mic_id, enz_id) in enumerate(zip(d.mic_ids, d.enzyme_ids)):
        er_ids = [
            er_id for er_id in S.columns
            if er_id.split(ID_SEPARATOR)[0] == enz_id
        ]
        out.loc[mic_id, er_ids] = i + 1
    return out.values.astype(int).tolist()


def get_conc_init(
    ms: MeasurementSet,
    kinetic_model: KineticModel,
    priors: PriorSet,
    config: MaudConfig
) -> StanData:
    """Get the initial mic concentrations for the ODE solver.

    The initial concentration for a measured mic is the measured value.

    The initial concentration for an unmeasured unbalanced mic in each
    experiment is its prior mean.

    The initial concentration for an unmeasured balanced mics is the default.

    :param mi: a MaudInput object

    """
    out = []
    y = ms.yconc.set_index("target_id")["value"]
    for mic in kinetic_model.mics:
        if mic.id in y.index:
            out.append(y.loc[mic.id])
        elif mic.balanced:
            out.append(config.default_initial_concentration)
        else:
            out.append(priors.conc_unbalanced.location.loc[mic.id])
    return StanData(out)


def unpack_priors(p: Union[IndPrior1d, IndPrior2d]) -> List[List[float]]:
    return [p.location.tolist(), p.scale.tolist()]
