import argparse
import os
from sbml_functions import read_sbml_file, StanReadySbmlModel
from typing import List


EXAMPLE_SBML_INPUT = '../data_in/t_brucei.xml'
PATH_TO_OUTPUT_DIRECTORY = '../stan/autogen'


def get_stan_program(m: StanReadySbmlModel) -> str:
    return '\n\n'.join([_build_functions(m),
                        _build_get_derived_quantities(m),
                        _build_get_fluxes_function(m),
                        _build_ode_function(m),
                        _build_steady_state_function()])

def _build_function_from_definition(name: str,
                                    arguments: List[str],
                                    body:str) -> str:
    arguments_string = ','.join([f'real {argument}' for argument in arguments])
    first_line = f"real {name}({arguments_string}){{"
    return_line = f"  return {body};"
    close_braces_line = "}"
    return '\n'.join([first_line, return_line, close_braces_line])
    

def _build_functions(m: StanReadySbmlModel) -> str:
    return '\n\n'.join([
        _build_function_from_definition(name,
                                        definition['arguments'],
                                        definition['body'])
        for name, definition in m.function_definitions.items()
    ])

def _build_get_derived_quantities(m: StanReadySbmlModel) -> str:
    first_line = "real[] get_derived_quantities(real[] ode_metabolites, real[] known_reals){"
    known_real_lines = [
        f"  real {kr} = known_reals[{i+1}];"
        for i, kr in enumerate(m.known_reals.keys())
    ]
    ode_metabolite_lines = [
        f"  real {metabolite} = ode_metabolites[{i+1}];"
        for i, metabolite in enumerate(m.ode_metabolites)
    ]
    derived_quantity_lines = [
        f"  real {quantity} = {expression};"
        for quantity, expression in m.assignment_expressions.items()
    ]
    return_line = f"  return {{{', '.join(m.assignment_expressions.keys())}}};"
    close_braces_line = "}"
    return '\n'.join([first_line,
                      *known_real_lines,
                      *ode_metabolite_lines,
                      *derived_quantity_lines,
                      return_line,
                      close_braces_line])


def _build_get_fluxes_function(m: StanReadySbmlModel) -> str:
    first_block = """real[] get_fluxes(real[] ode_metabolites,
                    real[] kinetic_parameters,
                    real[] known_reals){"""
    known_reals_lines = [
        f"  real {kr} = known_reals[{i+1}];"
        for i, kr in enumerate(m.known_reals.keys())
    ]
    ode_metabolite_lines = [
        f"  real {metabolite} = ode_metabolites[{i+1}];"
        for i, metabolite in enumerate(m.ode_metabolites.keys())
    ]
    kinetic_parameter_lines = [
        f"  real {parameter} = kinetic_parameters[{i+1}];"
        for i, parameter in enumerate(m.kinetic_parameters.keys())
    ]
    derived_quantity_calculation_line = (
        f"  real derived_quantities[{len(m.assignment_expressions.keys())}] = "
        "get_derived_quantities(ode_metabolites, known_reals);"
    )
    derived_quantity_unpack_lines = [
        f"  real {quantity} = derived_quantities[{i+1}];"
        for i, quantity in enumerate(m.assignment_expressions.keys())
    ]
    kinetic_expression_lines = [
        f"  real {parameter} = {expression};"
        for parameter, expression in m.kinetic_expressions.items()
    ]
    return_line = f"  return {{{', '.join(m.kinetic_expressions.keys())}}};"
    close_braces_line = "}"
    return '\n'.join([first_block,
                      *known_reals_lines,
                      *ode_metabolite_lines,
                      *kinetic_parameter_lines,
                      derived_quantity_calculation_line,
                      *derived_quantity_unpack_lines,
                      *kinetic_expression_lines,
                      return_line,
                      close_braces_line])


def _build_ode_function(m: StanReadySbmlModel) -> str:
    definition_line = "real[] get_odes(real[] fluxes){"
    unpack_lines = [
        f"  real {flux} = fluxes[{i+1}];"
        for i, flux in enumerate(m.kinetic_expressions.keys())
    ]
    return_line_open = "  return {"
    return_body_lines = [
        f"    {expression},  // {flux}"
        for flux, expression in m.ode_expressions.items()
    ]
    return_body_lines[-1] = return_body_lines[-1].replace(',', '')
    return_close_line = "  };"
    close_braces_line = "}"
    return '\n'.join([definition_line,
                      *unpack_lines,
                      return_line_open,
                      *return_body_lines,
                      return_close_line,
                      close_braces_line])


def _build_steady_state_function():
    return """real[] steady_state_equation(real t,
                             real[] ode_metabolites,
                             real[] kinetic_parameters,
                             real[] known_reals,
                             int[] known_ints){
    return get_odes(get_fluxes(ode_metabolites, kinetic_parameters, known_reals));\n}"""


if __name__ == '__main__': 
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file',
                        type=str,
                        default=EXAMPLE_SBML_INPUT,
                        help='Sbml file to convert')
    args = parser.parse_args()

    print('Parsing sbml input...')
    parsed_sbml = read_sbml_file(args.input_file)
    print('Converting parsed sbml to Stan...')
    stan_program = get_stan_program(parsed_sbml)

    output_filename = args.input_file.split('/')[-1].replace('.xml', '.stan')
    output_path = os.path.join(os.path.dirname(__file__), '../stan/autogen/' + output_filename)
    print(f'Writing Stan code to {output_path}...')

    with open(output_path, 'w') as f:
        f.write(stan_program)
        f.close()
    print('Finished!')
