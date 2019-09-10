from enzymekat._data_model import (
    EnzymeKatInput,
    Metabolite,
    Reaction,
    Modifier,
    Parameter,
    LSPrior,
    Experiment,
    Measurement,
    Compartment
)
import toml

def load_enzymekat_input_from_toml(filepath: str, id: str ='eki') -> EnzymeKatInput:
    """
    Load an EnzymeKatInput object from a suitable toml file
    
    :param filepath: path to a toml file
    :param id: id for the output object

    """
    eki = EnzymeKatInput(id)
    parsed_toml = toml.load(filepath)
    for c in parsed_toml['compartments']:
        cmp = Compartment(
            id=c['id'],
            name=c['name'],
            volume=c['volume']
        )
        eki.compartments.update({cmp.id: cmp})
    for m in parsed_toml['metabolites']:
        met = Metabolite(
            id=m['id'],
            name=m['name'],
            balanced=m['balanced'],
            compartment=m['compartment']
        )
        eki.metabolites.update({met.id: met})
    for r in parsed_toml['reactions']:
        rxn_params = {}
        for p in r['parameters']:
            param_met = p['metabolite'] if 'metabolite' in p.keys() else None
            param_prior = LSPrior(p['location'], p['scale'], 'lognormal')
            rxn_params.update(
                {p['id']: Parameter(p['id'], r['id'], param_met, param_prior)}
            )
        if 'modifiers' in r.keys():
            rxn_modifiers = {}
            for modifier_type, mets in r['modifiers'].items():
                for met in mets:
                    rxn_modifiers[met] = Modifier(
                        metabolite=met, modifier_type=modifier_type
                    )
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
        eki.reactions.update({rxn.id: rxn})
    for e in parsed_toml['experiments']:
        experiment = Experiment(id=e['id'])
        measurements = {}
        for measurement_type in ['metabolite', 'enzyme', 'reaction']:
            type_msmts = {}
            for m in e['measurements'][measurement_type]:
                msmt = Measurement(
                    value=m['value'],
                    uncertainty=m['uncertainty'],
                    scale='ln',
                    measurement_type=measurement_type
                )
                type_msmts.update({m['target']: msmt})
            measurements[measurement_type] = type_msmts
        experiment.met_meas.update(measurements['metabolite'])
        experiment.rxn_meas.update(measurements['reaction'])
        experiment.enz_meas.update(measurements['enzyme'])
        for ump in e['unbalanced_metabolite_priors']:
            prior = LSPrior(location=ump['location'], scale=ump['scale'], distribution='lognormal')
            experiment.unb_met_priors.update({ump['metabolite_id']: prior})
        eki.experiments.update({experiment.id: experiment})
    return eki

