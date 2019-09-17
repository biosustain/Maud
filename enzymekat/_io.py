from enzymekat._data_model import (
    EnzymeKatInput,
    KineticModel,
    Metabolite,
    Reaction,
    Modifier,
    Parameter,
    Prior,
    Experiment,
    Measurement,
    Compartment
)
import toml

MECHANISM_TO_PARAM_IDS = {
    'uniuni': ['Keq', 'Kcat1', 'Kcat2', 'Ka']
}

def load_enzymekat_input_from_toml(filepath: str, id: str ='eki') -> EnzymeKatInput:
    """
    Load an EnzymeKatInput object from a suitable toml file
    
    :param filepath: path to a toml file
    :param id: id for the output object

    """
    kinetic_model = KineticModel(id)
    parsed_toml = toml.load(filepath)
    for c in parsed_toml['compartments']:
        cmp = Compartment(
            id=c['id'],
            name=c['name'],
            volume=c['volume']
        )
        kinetic_model.compartments.update({cmp.id: cmp})
    for m in parsed_toml['metabolites']:
        met = Metabolite(
            id=m['id'],
            name=m['name'],
            balanced=m['balanced'],
            compartment=m['compartment']
        )
        kinetic_model.metabolites.update({met.id: met})
    for r in parsed_toml['reactions']:
        rxn_params = {}
        for param_id in MECHANISM_TO_PARAM_IDS[r['mechanism']]:
            rxn_params.update({param_id: Parameter(param_id, r['id'])})
        if 'allosteric_inhibitors' in r.keys():
            rxn_modifiers = {}
            for met in r['allosteric_inhibitors']:
                rxn_modifiers[met] = Modifier(met, 'allosteric_inhibition')
        else:
            rxn_modifiers = None
        rxn = Reaction(
            id=r['id'],
            name=r['name'] if r['name'] else None,
            reversible=r['reversible'] if 'reversible' in r.keys() else None,
            is_exchange=r['is_exchange'] if 'is_exchange' in r.keys() else None,
            stoichiometry=r['stoichiometry'],
            modifiers=rxn_modifiers,
            parameters=rxn_params,
            rate_law=r['rate_law'] if 'rate_law' in r.keys() else None,
            enzymes=r['enzymes'] if 'enzymes' in r.keys() else None 
        )
        kinetic_model.reactions.update({rxn.id: rxn})
    experiments = {}
    for e in parsed_toml['experiments']:
        experiment = Experiment(id=e['id'])
        measurements = {}
        for target_type in ['metabolite', 'enzyme', 'reaction']:
            type_msmts = {}
            for m in e[target_type + '_measurements']:
                msmt = Measurement(
                    target_id=m['target_id'],
                    value=m['value'],
                    uncertainty=m['uncertainty'],
                    scale='ln',
                    target_type=target_type
                )
                type_msmts.update({m['target_id']: msmt})
            measurements[target_type] = type_msmts
        experiment.met_meas.update(measurements['metabolite'])
        experiment.rxn_meas.update(measurements['reaction'])
        experiment.enz_meas.update(measurements['enzyme'])
        experiments.update({experiment.id: experiment})
    priors = {}
    for um_id, umps in parsed_toml['priors']['unbalanced_metabolites'].items():
        for ump in umps:
            prior_id = f"{um_id}_{ump['target_id']}"
            priors[prior_id] = Prior(
                id = prior_id,
                target_id = ump['target_id'],
                location=ump['location'],
                scale=ump['scale'],
                target_type='unbalanced_metabolite'
            )
    for kp_id, kpps in parsed_toml['priors']['kinetic_parameters'].items():
        for kpp in kpps:
            prior_id = f"{kp_id}_{kpp['target_id']}"
            priors[prior_id] = Prior(
                id = prior_id,
                target_id = kpp['target_id'],
                location=kpp['location'],
                scale=kpp['scale'],
                target_type='kinetic_parameter'
            )
    return EnzymeKatInput(kinetic_model=kinetic_model,
                          priors=priors,
                          experiments=experiments)

