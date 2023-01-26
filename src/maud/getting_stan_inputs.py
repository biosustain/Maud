"""Provides function get_stan_inputs for generating Stan input dictionaries."""

from dataclasses import fields
from typing import Dict, Iterable, List, Tuple, Union

from scipy.stats import gmean

from maud.data_model.experiment import Experiment, MeasurementType
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
    PhosphorylationModifyingEnzyme,
    Reaction,
    ReactionMechanism,
)
from maud.data_model.maud_config import MaudConfig
from maud.data_model.maud_parameter import ParameterSet
from maud.data_model.prior import IndPrior1d, IndPrior2d


def get_stan_inputs(
    parameters: ParameterSet,
    experiments: List[Experiment],
    kinetic_model: KineticModel,
    config: MaudConfig,
) -> Tuple[Dict, Dict]:
    """Get inputs to model.stan and prediction_model.stan.

    :param ms: a MeasurementSet object
    :param priors: a PriorSet object
    :param kinetic_model: a KineticModel object
    :param config: a MaudConfig object

    """
    network_properties_input = get_network_properties_input(
        kinetic_model, parameters
    )
    priors_input_train, priors_input_test = get_prior_inputs(parameters)
    experiments_input_train, experiments_input_test = get_experiments_input(
        experiments, kinetic_model
    )
    config_input = get_config_input(config)
    conc_init_train, conc_init_test = get_conc_init(
        experiments, kinetic_model, config
    )
    return (
        network_properties_input
        | priors_input_train
        | experiments_input_train
        | config_input
        | {"conc_init": conc_init_train},
        network_properties_input
        | priors_input_test
        | experiments_input_test
        | config_input
        | {"conc_init": conc_init_test},
    )


def get_prior_inputs(parameters: ParameterSet) -> Tuple[Dict, Dict]:
    """Get the priors component of an input to Maud's Stan model."""
    ind_priors_train = {
        f"priors_{p.name}": [p.prior.location, p.prior.scale]
        for p in map(lambda f: getattr(parameters, f.name), fields(parameters))
        if p.prior_in_train_model
        and (isinstance(p.prior, IndPrior1d) or isinstance(p.prior, IndPrior2d))
    }
    ind_priors_test = {
        f"priors_{p.name}": [p.prior.location, p.prior.scale]
        for p in map(lambda f: getattr(parameters, f.name), fields(parameters))
        if p.prior_in_test_model
        and (isinstance(p.prior, IndPrior1d) or isinstance(p.prior, IndPrior2d))
    }
    assert hasattr(parameters.dgf.prior, "covariance_matrix")
    dgf_priors = {
        "prior_loc_dgf": parameters.dgf.prior.location,
        "prior_cov_dgf": parameters.dgf.prior.covariance_matrix,
    }
    return (ind_priors_train | dgf_priors, ind_priors_test)


def get_network_properties_input(
    kinetic_model: KineticModel, parameters: ParameterSet
) -> Dict:
    """Get a dictionary of Stan data variables representing a kinetic model.

    :param kinetic_model: A KineticModel object
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
    km_ids = parameters.km.ids[0]
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
            prod_km_ids = (
                [ID_SEPARATOR.join([edge.enzyme_id, p]) for p in prod_ids]
                if rxn.mechanism
                != ReactionMechanism.irreversible_michaelis_menten
                else []
            )
        sub_code_by_edge.append([mic_codes[s] for s in sub_ids])
        sub_km_code_by_edge.append([km_codes[k] for k in sub_km_ids])
        prod_code_by_edge.append([mic_codes[p] for p in prod_ids])
        prod_km_code_by_edge.append([km_codes[k] for k in prod_km_ids])

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
        allosteric_enzyme_ids = list(
            set(a.enzyme_id for a in kinetic_model.allosteries)
        )
        allostery_type = [
            int(a.modification_type.value) for a in kinetic_model.allosteries
        ]
        allostery_mic = [mic_codes[a.mic_id] for a in kinetic_model.allosteries]
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
        allosteric_enzyme_ids = []
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
    ci_ix_long, ci_ix_bounds = encode_ragged(ci_by_edge)
    allostery_ix_long, allostery_ix_bounds = encode_ragged(allostery_by_edge)
    phosphorylation_ix_long, phosphorylation_ix_bounds = encode_ragged(
        phosphorylation_by_edge
    )
    return {
        "N_mic": len(kinetic_model.mics),
        "N_edge_sub": len(sub_by_edge_long),
        "N_edge_prod": len(prod_by_edge_long),
        "N_edge": len(S.columns),
        "N_unbalanced": len(unbalanced_mic_codes),
        "N_enzyme": len(enzymes),
        "N_phosphorylation": len(phosphorylation_codes),
        "N_pme": len(pme_codes),
        "N_competitive_inhibition": len(ci_ix_long),
        "N_allostery": len(allostery_ix_long),
        "N_allosteric_enzyme": len(allosteric_enzyme_ids),
        "N_drain": len(drain_codes),
        "N_km": len(km_codes),
        "N_sub_km": len(sub_km_ix_by_edge_long),
        "N_prod_km": len(prod_km_ix_by_edge_long),
        "S": S.values.tolist(),
        "N_reaction": len(reactions),
        "N_metabolite": len(metabolite_codes),
        "balanced_mic_ix": balanced_mic_codes,
        "unbalanced_mic_ix": unbalanced_mic_codes,
        "ci_mic_ix": ci_mic_codes,
        "edge_type": edge_mechanism_code,
        "edge_to_enzyme": edge_enzyme_code,
        "edge_to_tc": edge_tc_code,
        "edge_to_drain": edge_drain_code,
        "edge_to_reaction": edge_reaction_code,
        "water_stoichiometry": edge_water_stoichiometry,
        "transported_charge": edge_transported_charge,
        "mic_to_met": mic_met_code,
        "subunits": [e.subunits for e in kinetic_model.enzymes],
        "sub_by_edge_long": sub_by_edge_long,
        "sub_by_edge_bounds": sub_by_edge_bounds,
        "prod_by_edge_long": prod_by_edge_long,
        "prod_by_edge_bounds": prod_by_edge_bounds,
        "sub_km_ix_by_edge_long": sub_km_ix_by_edge_long,
        "sub_km_ix_by_edge_bounds": sub_km_ix_by_edge_bounds,
        "prod_km_ix_by_edge_long": prod_km_ix_by_edge_long,
        "prod_km_ix_by_edge_bounds": prod_km_ix_by_edge_bounds,
        "ci_ix_long": ci_ix_long,
        "ci_ix_bounds": ci_ix_bounds,
        "allostery_ix_long": allostery_ix_long,
        "allostery_ix_bounds": allostery_ix_bounds,
        "allostery_type": allostery_type,
        "allostery_mic": allostery_mic,
        "phosphorylation_ix_long": phosphorylation_ix_long,
        "phosphorylation_ix_bounds": phosphorylation_ix_bounds,
        "phosphorylation_type": phosphorylation_type,
        "phosphorylation_pme": phosphorylation_pme,
    }


def get_experiments_input(
    experiments: List[Experiment], kinetic_model: KineticModel
):
    """Get stan input data pertaining to experiments."""
    experiments_train = [e for e in experiments if e.is_train]
    experiments_test = [e for e in experiments if e.is_test]
    reaction_codes = codify_maud_object(kinetic_model.reactions)
    enzyme_codes = codify_maud_object(kinetic_model.enzymes)
    mic_codes = codify_maud_object(kinetic_model.mics)
    pme_codes = (
        codify_maud_object(kinetic_model.phosphorylation_modifying_enzymes)
        if kinetic_model.phosphorylation_modifying_enzymes is not None
        else dict()
    )
    experiment_codes_train = codify_maud_object(experiments_train)
    measurement_types = [
        MeasurementType.MIC,
        MeasurementType.ENZYME,
        MeasurementType.FLUX,
    ]
    m_conc_train, m_enz_train, m_flux_train = (
        [
            m
            for e in experiments
            for m in e.measurements
            if e.is_train and m.target_type == t
        ]
        for t in measurement_types
    )
    expids_m_conc_train, expids_m_enz_train, expids_m_flux_train = (
        [
            e.id
            for e in experiments
            for m in e.measurements
            if e.is_train and m.target_type == t
        ]
        for t in measurement_types
    )
    enz_ko_by_experiment_train, enz_ko_by_experiment_test = (
        [
            [enzyme_codes[eko.enzyme] for eko in experiment.enzyme_knockouts]
            if experiment.enzyme_knockouts is not None
            else []
            for experiment in exps
        ]
        for exps in (experiments_train, experiments_test)
    )
    pme_ko_by_experiment_train, pme_ko_by_experiment_test = (
        [
            [pme_codes[pko.pme] for pko in experiment.pme_knockouts]
            if experiment.pme_knockouts is not None
            else []
            for experiment in exps
        ]
        for exps in (experiments_train, experiments_test)
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
    input_train = {
        "N_experiment_train": len(experiments_train),
        "N_flux_measurement_train": len(m_flux_train),
        "N_enzyme_measurement_train": len(m_enz_train),
        "N_conc_measurement_train": len(m_conc_train),
        "N_enzyme_knockout_train": len(enz_knockout_long_train),
        "N_pme_knockout_train": len(pme_knockout_long_train),
        "temperature_train": [e.temperature for e in experiments_train],
        "enzyme_knockout_train_long": enz_knockout_long_train,
        "enzyme_knockout_train_bounds": enz_knockout_bounds_train,
        "pme_knockout_train_long": pme_knockout_long_train,
        "pme_knockout_train_bounds": pme_knockout_bounds_train,
        "yconc_train": [m.value for m in m_conc_train],
        "sigma_yconc_train": [m.error_scale for m in m_conc_train],
        "experiment_yconc_train": [
            experiment_codes_train[expid] for expid in expids_m_conc_train
        ],
        "mic_ix_yconc_train": [mic_codes[m.target_id] for m in m_conc_train],
        "yflux_train": [m.value for m in m_flux_train],
        "sigma_yflux_train": [m.error_scale for m in m_flux_train],
        "experiment_yflux_train": [
            experiment_codes_train[expid] for expid in expids_m_flux_train
        ],
        "reaction_yflux_train": [
            reaction_codes[m.target_id] for m in m_flux_train
        ],
        "yenz_train": [m.value for m in m_enz_train],
        "sigma_yenz_train": [m.error_scale for m in m_enz_train],
        "experiment_yenz_train": [
            experiment_codes_train[expid] for expid in expids_m_enz_train
        ],
        "enzyme_yenz_train": [enzyme_codes[m.target_id] for m in m_enz_train],
    }
    input_test = {
        "N_experiment_test": len(experiments_test),
        "N_enzyme_knockout_test": len(enz_knockout_long_test),
        "N_pme_knockout_test": len(pme_knockout_long_test),
        "temperature_test": [e.temperature for e in experiments_test],
        "enzyme_knockout_test_long": enz_knockout_long_test,
        "enzyme_knockout_test_bounds": enz_knockout_bounds_test,
        "pme_knockout_test_long": pme_knockout_long_test,
        "pme_knockout_test_bounds": pme_knockout_bounds_test,
    }
    return input_train, input_test


def get_config_input(config: MaudConfig):
    """Get Stan input related to algorithm configuration."""
    return {
        "likelihood": int(config.likelihood),
        "drain_small_conc_corrector": config.drain_small_conc_corrector,
        "reject_non_steady": int(config.reject_non_steady),
        "steady_state_threshold_abs": config.steady_state_threshold_abs,
        "steady_state_threshold_rel": config.steady_state_threshold_rel,
        "rel_tol": config.ode_config.rel_tol,
        "abs_tol": config.ode_config.abs_tol,
        "timepoint": config.ode_config.timepoint,
        "max_num_steps": config.ode_config.max_num_steps,
    }


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
    experiments: List[Experiment],
    kinetic_model: KineticModel,
    config: MaudConfig,
) -> Tuple[List[List[float]], List[List[float]]]:
    """Get the initial balanced mic concentrations for the ODE solver.

    The initial concentration for a measured mic is the average measured value.

    The initial concentration for an unmeasured mic is the default.

    """
    default = config.default_initial_concentration
    balanced_mics = [m for m in kinetic_model.mics if m.balanced]
    inits = [[default for _ in balanced_mics] for e in experiments]
    for i, experiment in enumerate(experiments):
        msmts = [
            m
            for m in experiment.measurements
            if m.target_type == MeasurementType.MIC
        ]
        for j, mic in enumerate(balanced_mics):
            if any(mic.id == m.target_id for m in msmts):
                inits[i][j] = gmean(
                    [m.value for m in msmts if mic.id == m.target_id]
                )
    inits_train = [
        inits[i] for i, exp in enumerate(experiments) if exp.is_train
    ]
    inits_test = [inits[i] for i, exp in enumerate(experiments) if exp.is_test]
    return inits_train, inits_test


CodifiableMaudObject = Union[
    Allostery,
    MetaboliteInCompartment,
    Metabolite,
    EnzymeReaction,
    Enzyme,
    Phosphorylation,
    PhosphorylationModifyingEnzyme,
    Reaction,
    Experiment,
    CompetitiveInhibition,
]


def codify_maud_object(lx: Iterable[CodifiableMaudObject]) -> Dict[str, int]:
    """Turn an iterable of Maud objects into a dictionary of Stan codes."""
    return dict(zip([x.id for x in lx], range(1, len(list(lx)) + 1)))
