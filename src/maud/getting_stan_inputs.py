"""Code for creating StanInput objects."""

from typing import Dict, Iterable, List, Optional, Tuple, Union

import numpy as np

from maud.data_model.hardcoding import ID_SEPARATOR
from maud.data_model.kinetic_model import (
    Allostery,
    CompetitiveInhibition,
    Edge,
    Enzyme,
    EnzymeReaction,
    KineticModel,
    Metabolite,
    MetaboliteInCompartment,
    Phosphorylation,
    PhosphorylationModifyingEnzyme,
    Reaction,
)
from maud.data_model.maud_config import MaudConfig
from maud.data_model.measurement_set import Experiment, MeasurementSet
from maud.data_model.prior_set import IndPrior1d, IndPrior2d, PriorSet
from maud.data_model.stan_input import StanInput


CodifiableMaudObject = Union[
    Allostery,
    MetaboliteInCompartment,
    Metabolite,
    EnzymeReaction,
    Enzyme,
    Edge,
    Phosphorylation,
    PhosphorylationModifyingEnzyme,
    Reaction,
    Experiment,
    CompetitiveInhibition,
]


def get_stan_inputs(
    ms: MeasurementSet,
    priors: PriorSet,
    kinetic_model: KineticModel,
    config: MaudConfig,
) -> Tuple[StanInput, Optional[StanInput]]:
    """Get inputs to model.stan and prediction_model.stan.

    :param ms: a MeasurementSet object
    :param priors: a PriorSet object
    :param kinetic_model: a KineticModel object
    :param config: a MaudConfig object

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
    experiments_train = [e for e in ms.experiments if e.is_train]
    experiments_test = [e for e in ms.experiments if e.is_test]
    experiment_ids_train = [exp.id for exp in experiments_train]
    experiment_ids_test = [exp.id for exp in experiments_test]
    experiment_codes_train = codify_maud_object(experiments_train)
    km_ids = priors.km.stan_variable.ids[0]
    km_codes = dict(zip(km_ids, range(1, len(km_ids) + 1)))
    if allosteries is not None:
        tc_codes = codify_maud_object(
            [
                e
                for e in enzymes
                if any(a.enzyme_id == e.id for a in allosteries)
            ]
        )
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
            kinetic_model.phosphorylations
        )
        assert kinetic_model.phosphorylation_modifying_enzymes is not None
        pme_codes = codify_maud_object(
            kinetic_model.phosphorylation_modifying_enzymes
        )
    else:
        phosphorylation_codes = dict()
        pme_codes = dict()
    drain_codes = codify_maud_object(kinetic_model.drains)
    edge_codes = codify_maud_object(edges)
    er_codes = codify_maud_object(kinetic_model.ers)
    reaction_codes = codify_maud_object(reactions)
    reaction_to_mechanism = {rxn.id: rxn.mechanism for rxn in reactions}
    edge_mechanism_code = [
        int(reaction_to_mechanism[e.id].value)
        if isinstance(e, Reaction)
        else int(reaction_to_mechanism[e.reaction_id].value)
        for e in edges
    ]
    edge_enzyme_code = [
        enzyme_codes[er.enzyme_id] if isinstance(er, EnzymeReaction) else 0
        for er in edges
    ]
    edge_er_code = [
        er_codes[er.id] if isinstance(er, EnzymeReaction) else 0
        for er in edges
    ]
    edge_drain_code = [
        drain_codes[d.id] if isinstance(d, Reaction) else 0 for d in edges
    ]
    edge_reaction_code = [
        reaction_codes[e.id]
        if isinstance(e, Reaction)
        else reaction_codes[e.reaction_id]
        for e in edges
    ]
    edge_water_stoichiometry = [
        reaction_by_id[e.reaction_id].water_stoichiometry
        if isinstance(e, EnzymeReaction)
        else 0
        for e in edges
    ]
    edge_transported_charge = [
        reaction_by_id[e.reaction_id].transported_charge
        if isinstance(e, EnzymeReaction)
        else 0
        for e in edges
    ]
    mic_met_code = [metabolite_codes[m.metabolite_id] for m in mics]
    balanced_mic_codes = [mic_codes[mic.id] for mic in mics if mic.balanced]
    unbalanced_mic_codes = [
        mic_codes[mic.id] for mic in mics if not mic.balanced
    ]
    ci_mic_codes = (
        [mic_codes[ci.mic_id] for ci in kinetic_model.competitive_inhibitions]
        if kinetic_model.competitive_inhibitions is not None
        else []
    )
    edge_tc_code = [
        tc_codes[e.enzyme_id]
        if isinstance(e, EnzymeReaction) and e.enzyme_id in tc_codes.keys()
        else 0
        for e in edges
    ]
    # ragged arrays
    sub_code_by_edge = []
    sub_km_code_by_edge = []
    prod_code_by_edge = []
    prod_km_code_by_edge = []
    for edge in edges:
        if isinstance(edge, Reaction):  # drain case
            rxn = edge
        else:
            rxn = next(r for r in reactions if r.id == edge.reaction_id)
        sub_ids = [
            s for s in rxn.stoichiometry.keys() if rxn.stoichiometry[s] < 0
        ]
        prod_ids = [
            p for p in rxn.stoichiometry.keys() if rxn.stoichiometry[p] > 0
        ]
        if isinstance(edge, Reaction):  # drain case
            sub_km_ids = []
            prod_km_ids = []
        else:
            sub_km_ids = [
                ID_SEPARATOR.join([edge.enzyme_id, s]) for s in sub_ids
            ]
            prod_km_ids = [
                ID_SEPARATOR.join([edge.enzyme_id, p]) for p in prod_ids
            ]
        sub_code_by_edge.append([mic_codes[s] for s in sub_ids])
        sub_km_code_by_edge.append([km_codes[k] for k in sub_km_ids])
        prod_code_by_edge.append([mic_codes[p] for p in prod_ids])
        prod_km_code_by_edge.append([km_codes[k] for k in prod_km_ids])
    enz_ko_by_experiment_train, enz_ko_by_experiment_test = (
        [
            [
                enzyme_codes[eko.enzyme_id]
                for eko in ms.enzyme_knockouts
                if eko.experiment_id == e
            ]
            for e in exps
        ]
        if ms.enzyme_knockouts is not None
        else [[] for e in exps]
        for exps in (experiments_train, experiments_test)
    )
    pme_ko_by_experiment_train, pme_ko_by_experiment_test = (
        [
            [
                pme_codes[pko.pme_id]
                for pko in ms.pme_knockouts
                if pko.experiment_id == e
            ]
            for e in exps
        ]
        if ms.pme_knockouts is not None
        else [[] for e in exps]
        for exps in (experiments_train, experiments_test)
    )
    ci_by_edge: List[List[int]] = [[] for er in kinetic_model.edges]
    if kinetic_model.competitive_inhibitions is not None:
        for ci in kinetic_model.competitive_inhibitions:
            edge_pos = next(
                i
                for i, e in enumerate(kinetic_model.edges)
                if (
                    isinstance(e, EnzymeReaction)
                    and e.enzyme_id == ci.enzyme_id
                )
            )
            ci_by_edge[edge_pos].append(ci_codes[ci.id])
    allostery_by_edge: List[List[int]] = [[] for er in kinetic_model.edges]
    if kinetic_model.allosteries is not None:
        allostery_type = [
            int(a.modification_type.value) for a in kinetic_model.allosteries
        ]
        allostery_mic = [
            mic_codes[a.mic_id] for a in kinetic_model.allosteries
        ]
        for a in kinetic_model.allosteries:
            edge_pos = next(  # is this right???
                i
                for i, er in enumerate(kinetic_model.edges)
                if (
                    isinstance(er, EnzymeReaction)
                    and er.enzyme_id == a.enzyme_id
                )
            )
            allostery_by_edge[edge_pos].append(allostery_codes[a.id])
    else:
        allostery_type = []
        allostery_mic = []
    phosphorylation_by_edge: List[List[int]] = [
        [] for er in kinetic_model.edges
    ]
    if kinetic_model.phosphorylations is not None:
        phosphorylation_type = [
            int(p.modification_type.value)
            for p in kinetic_model.phosphorylations
        ]
        phosphorylation_pme = [
            pme_codes[p.modifying_enzyme_id]
            for p in kinetic_model.phosphorylations
        ]
        for p in kinetic_model.phosphorylations:
            edge_pos = next(
                i
                for i, er in enumerate(kinetic_model.edges)
                if (
                    isinstance(er, EnzymeReaction)
                    and er.enzyme_id == p.modified_enzyme_id
                )
            )
            phosphorylation_by_edge[edge_pos].append(
                phosphorylation_codes[p.id]
            )
    else:
        phosphorylation_type = []
        phosphorylation_pme = []
    sub_by_edge_long, sub_by_edge_bounds = encode_ragged(sub_code_by_edge)
    prod_by_edge_long, prod_by_edge_bounds = encode_ragged(prod_code_by_edge)
    sub_km_ix_by_edge_long, sub_km_ix_by_edge_bounds = encode_ragged(
        sub_km_code_by_edge
    )
    prod_km_ix_by_edge_long, prod_km_ix_by_edge_bounds = encode_ragged(
        prod_km_code_by_edge
    )
    enz_knockout_long_train, enz_knockout_bounds_train = encode_ragged(
        enz_ko_by_experiment_train
    )
    enz_knockout_long_test, enz_knockout_bounds_test = encode_ragged(
        enz_ko_by_experiment_test
    )
    pme_knockout_long_train, pme_knockout_bounds_train = encode_ragged(
        pme_ko_by_experiment_train
    )
    pme_knockout_long_test, pme_knockout_bounds_test = encode_ragged(
        pme_ko_by_experiment_test
    )
    ci_ix_long, ci_ix_bounds = encode_ragged(ci_by_edge)
    allostery_ix_long, allostery_ix_bounds = encode_ragged(allostery_by_edge)
    phosphorylation_ix_long, phosphorylation_ix_bounds = encode_ragged(
        phosphorylation_by_edge
    )
    promiscuous_enzymes = [
        e
        for e in kinetic_model.enzymes
        if len([er for er in kinetic_model.ers if er.enzyme_id == e.id]) > 1
    ]
    edge_by_promiscuous_enzyme = [
        [
            edge_codes[ID_SEPARATOR.join([er.enzyme_id, er.reaction_id])]
            for er in kinetic_model.ers
            if er.enzyme_id == pe.id
        ]
        for pe in promiscuous_enzymes
    ]
    (
        edge_by_promiscuous_enzyme_long,
        edge_by_promiscuous_enzyme_bounds,
    ) = encode_ragged(edge_by_promiscuous_enzyme)
    yconc_train, yenz_train, yflux_train = (
        df.loc[
            lambda df: df["experiment_id"].isin(
                [e.id for e in experiments_train]
            )
        ].reset_index()
        for df in (ms.yconc, ms.yenz, ms.yflux)
    )
    priors_conc_pme_train, priors_conc_pme_test = [
        unpack_priors_2d(priors.conc_pme, exp_ids=exp_ids)
        for exp_ids in (experiment_ids_train, experiment_ids_test)
    ]
    priors_conc_unbalanced_train, priors_conc_unbalanced_test = [
        unpack_priors_2d(priors.conc_unbalanced, exp_ids=exp_ids)
        for exp_ids in (experiment_ids_train, experiment_ids_test)
    ]
    priors_conc_enzyme_train, priors_conc_enzyme_test = [
        unpack_priors_2d(priors.conc_enzyme, exp_ids=exp_ids)
        for exp_ids in (experiment_ids_train, experiment_ids_test)
    ]
    priors_drain_train, priors_drain_test = [
        unpack_priors_2d(priors.drain, exp_ids=exp_ids)
        for exp_ids in (experiment_ids_train, experiment_ids_test)
    ]
    priors_pepism_train, priors_pepism_test = [
        unpack_priors_2d(
            priors.promiscuous_enzyme_proportion_invsm, exp_ids=exp_ids
        )
        for exp_ids in (experiment_ids_train, experiment_ids_test)
    ]
    stan_input_train = StanInput(
        dict(
            # dimensions
            N_mic=S.shape[0],
            N_edge=S.shape[1],
            N_reaction=len(reactions),
            N_enzyme=len(enzymes),
            N_enzyme_reaction=len(er_codes),
            N_edge_sub=len(sub_by_edge_long),
            N_edge_prod=len(prod_by_edge_long),
            N_unbalanced=len(unbalanced_mic_codes),
            N_metabolite=len(metabolite_codes),
            N_km=len(km_ids),
            N_sub_km=len(sub_km_ix_by_edge_long),
            N_prod_km=len(prod_km_ix_by_edge_long),
            N_drain=len(drain_codes),
            N_allostery=len(allosteries) if allosteries is not None else 0,
            N_allosteric_enzyme=len(tc_codes),
            N_phosphorylation=len(phosphorylation_codes),
            N_pme=len(pme_codes),
            N_competitive_inhibition=len(priors.ki.location),
            N_experiment=len(experiments_train),
            N_flux_measurement=len(yflux_train),
            N_enzyme_measurement=len(yenz_train),
            N_conc_measurement=len(yconc_train),
            N_enzyme_knockout=len(enz_knockout_long_train),
            N_pme_knockout=len(pme_knockout_long_train),
            N_promiscuous_enzyme=len(edge_by_promiscuous_enzyme),
            N_promiscuous_enzyme_edge=len(edge_by_promiscuous_enzyme_long),
            # measurements
            experiment_yconc=yconc_train["experiment_id"]
            .map(experiment_codes_train)
            .tolist(),
            mic_ix_yconc=yconc_train["target_id"].map(mic_codes).tolist(),
            yconc=yconc_train["measurement"].tolist(),
            sigma_conc=yconc_train["error_scale"].tolist(),
            experiment_yflux=yflux_train["experiment_id"]
            .map(experiment_codes_train)
            .to_list(),
            reaction_yflux=yflux_train["target_id"]
            .map(reaction_codes)
            .tolist(),
            yflux=yflux_train["measurement"].tolist(),
            sigma_flux=yflux_train["error_scale"].tolist(),
            experiment_yenz=yenz_train["experiment_id"]
            .map(experiment_codes_train)
            .tolist(),
            enzyme_yenz=yenz_train["target_id"].map(enzyme_codes).tolist(),
            yenz=yenz_train["measurement"].tolist(),
            sigma_enz=yenz_train["error_scale"].tolist(),
            # hardcoded priors
            prior_loc_dgf=priors.dgf.location.to_list(),
            prior_cov_dgf=priors.dgf.covariance_matrix.values.tolist(),
            priors_kcat=unpack_priors_1d(priors.kcat),
            priors_km=unpack_priors_1d(priors.km),
            priors_ki=unpack_priors_1d(priors.ki),
            priors_dissociation_constant=unpack_priors_1d(
                priors.dissociation_constant
            ),
            priors_transfer_constant=unpack_priors_1d(
                priors.transfer_constant
            ),
            priors_kcat_pme=unpack_priors_1d(priors.kcat_pme),
            priors_psi=unpack_priors_1d(priors.psi, ids=experiment_ids_train),
            priors_conc_pme=priors_conc_pme_train,
            priors_conc_unbalanced=priors_conc_unbalanced_train,
            priors_conc_enzyme=priors_conc_enzyme_train,
            priors_drain=priors_drain_train,
            priors_promiscuous_enzyme_proportion_invsm=priors_pepism_train,
            # network properties
            S=S.values.tolist(),
            balanced_mic_ix=balanced_mic_codes,
            unbalanced_mic_ix=unbalanced_mic_codes,
            ci_mic_ix=ci_mic_codes,
            edge_type=edge_mechanism_code,
            edge_to_enzyme=edge_enzyme_code,
            edge_to_enzyme_reaction=edge_er_code,
            edge_to_tc=edge_tc_code,
            edge_to_drain=edge_drain_code,
            edge_to_reaction=edge_reaction_code,
            allostery_type=allostery_type,
            allostery_mic=allostery_mic,
            phosphorylation_type=phosphorylation_type,
            phosphorylation_pme=phosphorylation_pme,
            edge_by_promiscuous_enzyme_long=edge_by_promiscuous_enzyme_long,
            edge_by_promiscuous_enzyme_bounds=edge_by_promiscuous_enzyme_bounds,
            water_stoichiometry=edge_water_stoichiometry,
            transported_charge=edge_transported_charge,
            mic_to_met=mic_met_code,
            subunits=[e.subunits for e in kinetic_model.enzymes],
            temperature=[e.temperature for e in experiments_train],
            sub_by_edge_long=sub_by_edge_long,
            sub_by_edge_bounds=sub_by_edge_bounds,
            prod_by_edge_long=prod_by_edge_long,
            prod_by_edge_bounds=prod_by_edge_bounds,
            sub_km_ix_by_edge_long=sub_km_ix_by_edge_long,
            sub_km_ix_by_edge_bounds=sub_km_ix_by_edge_bounds,
            prod_km_ix_by_edge_long=prod_km_ix_by_edge_long,
            prod_km_ix_by_edge_bounds=prod_km_ix_by_edge_bounds,
            ci_ix_long=ci_ix_long,
            ci_ix_bounds=ci_ix_bounds,
            allostery_ix_long=allostery_ix_long,
            allostery_ix_bounds=allostery_ix_bounds,
            phosphorylation_ix_long=phosphorylation_ix_long,
            phosphorylation_ix_bounds=phosphorylation_ix_bounds,
            enzyme_knockout_long=enz_knockout_long_train,
            enzyme_knockout_bounds=enz_knockout_bounds_train,
            pme_knockout_long=pme_knockout_long_train,
            pme_knockout_bounds=pme_knockout_bounds_train,
            # configuration
            conc_init=get_conc_init(
                ms, kinetic_model, priors, config, mode="train"
            ),
            likelihood=int(config.likelihood),
            drain_small_conc_corrector=int(config.drain_small_conc_corrector),
            reject_non_steady=int(config.reject_non_steady),
            steady_state_threshold_abs=config.steady_state_threshold_abs,
            steady_state_threshold_rel=config.steady_state_threshold_rel,
            rel_tol=config.ode_config.rel_tol,
            abs_tol=config.ode_config.abs_tol,
            timepoint=config.ode_config.timepoint,
            max_num_steps=config.ode_config.max_num_steps,
        )
    )
    stan_input_test = (
        None
        if len(experiments_test) == 0
        else StanInput(
            dict(
                # priors
                priors_conc_pme=priors_conc_pme_test,
                priors_conc_unbalanced=priors_conc_unbalanced_test,
                priors_conc_enzyme=priors_conc_enzyme_test,
                priors_drain=priors_drain_test,
                # network properties
                S=S.values.tolist(),
                N_metabolite=len(kinetic_model.metabolites),
                N_km=len(km_ids),
                N_reaction=len(reactions),
                N_allosteric_enzyme=len(tc_codes),
                balanced_mic_ix=balanced_mic_codes,
                unbalanced_mic_ix=unbalanced_mic_codes,
                ci_mic_ix=ci_mic_codes,
                edge_type=edge_mechanism_code,
                edge_to_enzyme=edge_enzyme_code,
                edge_to_tc=edge_tc_code,
                edge_to_drain=edge_drain_code,
                edge_to_reaction=edge_reaction_code,
                water_stoichiometry=edge_water_stoichiometry,
                mic_to_met=mic_met_code,
                subunits=[e.subunits for e in kinetic_model.enzymes],
                temperature=[e.temperature for e in experiments_test],
                sub_by_edge_long=sub_by_edge_long,
                sub_by_edge_bounds=sub_by_edge_bounds,
                prod_by_edge_long=prod_by_edge_long,
                prod_by_edge_bounds=prod_by_edge_bounds,
                sub_km_ix_by_edge_long=sub_km_ix_by_edge_long,
                sub_km_ix_by_edge_bounds=sub_km_ix_by_edge_bounds,
                prod_km_ix_by_edge_long=prod_km_ix_by_edge_long,
                prod_km_ix_by_edge_bounds=prod_km_ix_by_edge_bounds,
                ci_ix_long=ci_ix_long,
                ci_ix_bounds=ci_ix_bounds,
                allostery_ix_long=allostery_ix_long,
                allostery_ix_bounds=allostery_ix_bounds,
                allostery_type=allostery_type,
                allostery_mic=allostery_mic,
                phosphorylation_ix_long=phosphorylation_ix_long,
                phosphorylation_ix_bounds=phosphorylation_ix_bounds,
                phosphorylation_type=phosphorylation_type,
                phosphorylation_pme=phosphorylation_pme,
                enzyme_knockout_long=enz_knockout_long_test,
                enzyme_knockout_bounds=enz_knockout_bounds_test,
                pme_knockout_long=pme_knockout_long_test,
                pme_knockout_bounds=pme_knockout_bounds_test,
                # configuration
                conc_init=get_conc_init(
                    ms, kinetic_model, priors, config, mode="test"
                ),
                likelihood=int(config.likelihood),
                drain_small_conc_corrector=int(
                    config.drain_small_conc_corrector
                ),
                rel_tol=config.ode_config.rel_tol,
                abs_tol=config.ode_config.abs_tol,
                timepoint=config.ode_config.timepoint,
                max_num_steps=config.ode_config.max_num_steps,
            )
        )
    )
    return stan_input_train, stan_input_test


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


def get_conc_init(
    ms: MeasurementSet,
    kinetic_model: KineticModel,
    priors: PriorSet,
    config: MaudConfig,
    mode: str,
) -> List[List[float]]:
    """Get the initial mic concentrations for the ODE solver.

    The initial concentration for a measured mic is the measured value.

    The initial concentration for an unmeasured unbalanced mic in each
    experiment is its prior mean.

    The initial concentration for an unmeasured balanced mics is the default.

    :param mi: a MaudInput object

    :param mode: either "train" or "test"

    """
    if mode not in ["train", "test"]:
        raise ValueError("Param mode must be either 'train' or 'test'.")
    out = []
    y = ms.yconc.set_index(["target_id", "experiment_id"])["measurement"]
    for exp in ms.experiments:
        if (mode == "train" and exp.is_train) or (
            mode == "test" and exp.is_test
        ):
            init_exp = []
            for mic in kinetic_model.mics:
                if (mic.id, exp.id) in y.index:
                    init_exp.append(y.loc[(mic.id, exp.id)])
                elif mic.balanced:
                    init_exp.append(config.default_initial_concentration)
                else:
                    init_exp.append(
                        np.exp(
                            priors.conc_unbalanced.location.loc[exp.id, mic.id]
                        )
                    )
            out.append(init_exp)
    return out


def unpack_priors_2d(
    p: IndPrior2d, exp_ids: Optional[List[str]]
) -> List[List[List[float]]]:
    """Turn an IndPrior2d object into a json-compatible list."""
    if exp_ids is None:
        exp_ids = p.location.index.tolist()
    if len(exp_ids) == 0:
        exp_ids = p.location.index.tolist()
    assert isinstance(exp_ids, List)
    loc = p.location.loc[exp_ids].values.tolist()
    scale = p.scale.loc[exp_ids].values.tolist()
    if len(loc) == 0:
        loc = [[] for _ in range(len(exp_ids))]
        scale = [[] for _ in range(len(exp_ids))]
    return [loc, scale]


def unpack_priors_1d(p: IndPrior1d, ids=None) -> List[List[float]]:
    """Turn an IndPrior1d object into a json-compatible list."""
    if ids is None:
        return [p.location.to_list(), p.scale.to_list()]
    else:
        return [p.location.loc[ids].to_list(), p.scale.loc[ids].to_list()]


def codify_maud_object(lx: Iterable[CodifiableMaudObject]) -> Dict[str, int]:
    """Turn an iterable of Maud objects into a dictionary of Stan codes."""
    return dict(zip([x.id for x in lx], range(1, len(list(lx)) + 1)))
