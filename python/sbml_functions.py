from collections import namedtuple
import pandas as pd
import libsbml
from libsbml import formulaToString
from pprint import pprint

StanReadySbmlModel = namedtuple('OdeSystem', ['ode_metabolites',
                                              'known_reals',
                                              'kinetic_parameters',
                                              'assignment_rules',
                                              'function_definitions',
                                              'kinetic_expressions',
                                              'ode_expressions'])


def read_sbml_file(f):
    reader = libsbml.SBMLReader()
    d = reader.readSBML(f)
    m = d.getModel()
    return StanReadySbmlModel(_get_ode_metabolites(m),
                              _get_known_reals(m),
                              _get_kinetic_parameters(m),
                              _get_assignment_rules(m),
                              _get_function_definitions(m),
                              _get_kinetic_expressions(m),
                              _get_ode_expressions(m))


def _get_ode_metabolites(m):
    assignment_rules = _get_assignment_rules(m)
    return {s.getId(): s.getInitialConcentration()
            for s in m.getListOfSpecies()
            if s.getId() not in assignment_rules.keys()}
    

def _get_known_reals(m):
    compartments = m.getListOfCompartments()
    compartment_volumes = {compartment.getId(): compartment.getVolume()
                           for compartment in compartments}
    fixed_parameters = {parameter.getId(): parameter.getValue()
                        for parameter in m.getListOfParameters()
                        if parameter.getConstant()}
    return {**compartment_volumes,
            **fixed_parameters}


def _get_kinetic_parameters(m):
    out = {}
    reactions = m.getListOfReactions()
    for reaction in reactions:
        parameters = reaction.getKineticLaw().getListOfParameters()
        for parameter in parameters:
            out[parameter.getId()] = parameter.getValue()
    return out


def _get_assignment_rules(m):
    rules = filter(lambda r: r.isAssignment(), m.getListOfRules())
    return {rule.getId(): rule.getFormula()
            for rule in rules}


def _get_function_definitions(m):
    out = {}
    for definition in m.getListOfFunctionDefinitions():
        arguments = [definition.getArgument(ix).getName()
                     for ix in range(definition.getNumArguments())]
        body = formulaToString(definition.getBody())
        out[definition.getId()] = {'arguments': arguments, 'body': body}
    return out
        
        

def _get_kinetic_expressions(m):
    return {reaction.getId(): reaction.getKineticLaw().getFormula()
            for reaction in m.getListOfReactions()}


def _get_ode_expressions(m):
    """Inspired by this code:
       https://www.dropbox.com/s/2bfpiausejp0gd0/convert_reactions.py?dl=0
    """

    def update_expressions(list_of_expressions, reactant_ref, species, reaction):
        out = list_of_expressions.copy()
        if reactant_ref.getSpecies() == species.getId():
            stoichiometry = reactant_ref.getStoichiometry()
            reaction_id = reaction.getId()
            out.append(f'{stoichiometry}*{reaction_id}')
        return out

    out = {}
    specieses = m.getListOfSpecies()
    reactions = m.getListOfReactions()
    ode_metabolites = list(_get_ode_metabolites(m).keys())
    for species in specieses:
        if species.getId() not in _get_assignment_rules(m).keys():
            reactant_expressions = []
            product_expressions = []
            for reaction in reactions:
                for reactant_ref_ix in range(reaction.getNumReactants()):
                    reactant_ref = reaction.getReactant(reactant_ref_ix)
                    reactant_expressions = update_expressions(reactant_expressions,
                                                            reactant_ref,
                                                            species,
                                                            reaction)
                for product_ref_ix in range(reaction.getNumProducts()):
                    product_ref = reaction.getProduct(product_ref_ix)
                    product_expressions = update_expressions(product_expressions,
                                                            product_ref,
                                                            species,
                                                            reaction)
            expr = ('+'.join(product_expressions)
                    + ''.join([f'-{i}' for i in reactant_expressions]))
            if expr != '':
                out[species.getId()] = expr
    return out
