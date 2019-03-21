import re

INPUT_FILE = "Ecoli glycolisis.mmd"
OUTPUT_FILE = "data/steady_state_autogen.stan"


def remove_non_numeric_bits(s, return_type=float):
    """Remove non-numeric bits from a string"""
    regex = "[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?"
    finds = re.findall(regex, s)
    if return_type == float:
        return float(finds[0])
    elif return_type == str:
        return finds[0]
    elif return_type == list:
        return finds
    else:
        raise ValueError('unsupported return type: ' + str(return_type))
    
def get_all_metabolites(input_file):
    with open(input_file, 'r') as f:
        for l in f:
            if 'concentration of metabolite' in l:
                yield re.findall("concentration of metabolite '(.*?)'", l)[0]


def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def get_named_quantity(input_file, regex, remove_non_numeric=True):
    out = dict()
    with open(input_file, 'r') as f:
        for l in f:
            if re.match('.*' + regex + '.*', l) is not None:
                split_line = l.split('=')
                name = split_line[0].strip()
                quantity = split_line[1].split(';')[0].replace('\t', '').strip()
                # turn square bracket indexes into subscripts
                name = re.sub(r"\[([^{}]+)\]", r"_\1", name)
                quantity = re.sub(r"\[([^{}]+)\]", r"_\1", quantity)
                if remove_non_numeric:
                    quantity = remove_non_numeric_bits(quantity)
                out[name] = quantity
    return out


def get_ode_metabolites(input_file):
    return {
        k.replace('init ', ''): float(v)
        for k, v in get_named_quantity(input_file, '; metabolite.*reactions', False).items()
        if is_float(v)
    }


def get_derived_quantities(input_file):
    derived_metabolites = {
        k: v
        for k, v in get_named_quantity(input_file, '; metabolite', False).items()
        if not is_float(v)
    }
    model_entities = get_named_quantity(input_file, 'model entity', False)
    concentrations_reaction = get_named_quantity(
        input_file, 'concentration of metabolite .*reactions', False
    )
    concentrations_fixed = get_named_quantity(
        input_file, 'concentration of metabolite .*fixed', False
    )
    concentrations_assignment = get_named_quantity(
        input_file, 'concentration of metabolite .*assignment', False

    )
    
    return {**derived_metabolites,
            **concentrations_reaction,
            **concentrations_fixed,
            **model_entities,
            **concentrations_assignment}
        

def get_known_reals(input_file):
    conserved_totals = get_named_quantity(INPUT_FILE, 'conserved total')
    compartments = get_named_quantity(INPUT_FILE, 'compartment')
    fixed_global = get_named_quantity(INPUT_FILE, 'global quantity .*fixed')
    fixed_compartment = get_named_quantity(INPUT_FILE, 'compartment .*fixed')
    fixed_metabolite = get_named_quantity(INPUT_FILE, '; metabolite .*fixed')
    return {**conserved_totals,
            **compartments,
            **fixed_global,
            **fixed_compartment,
            **fixed_metabolite}


def get_kinetic_functions(input_file):
    out = dict()
    with open(input_file, 'r') as f:
        for l in f:
            if 'kinetic function' in l:
                reaction = l.split(' ')[0].replace('FunctionFor', '')
                expression = l.split(' ')[2].replace('\t', '')
                out[reaction] = expression
    return out


def get_odes(input_file):
    out = dict()
    with open(input_file, 'r') as f:
        for l in f:
            if l[:4] == 'd/dt':
                reaction = re.findall('d/dt\((.*?)\)', l)[0]
                expression = (l.split(' ')[2]
                              .replace('\t', '')
                              .replace('FunctionFor', '')
                              .replace(';', ',')
                              .replace('*cell_cytoplasm', ''))  #TODO: deal with compartments properly
                out[reaction] = expression
    return out


def build_derived_quantity_function(known_reals,
                                    ode_metabolites,
                                    derived_quantity_expressions):
    first_line = """real[] get_derived_quantities(vector ode_metabolites, real[] known_reals){"""
    known_reals_lines = [
        f"  real {kr} = known_reals[{i+1}];"
        for i, kr in enumerate(known_reals.keys())
    ]
    ode_metabolite_lines = [
        f"  real {metabolite} = ode_metabolites[{i+1}];"
        for i, metabolite in enumerate(ode_metabolites)
    ]
    derivation_lines = [
        f"  real {quantity} = {expression};"
        for quantity, expression in derived_quantity_expressions.items()
    ]
    return_line = f"  return {{{', '.join(derived_quantity_expressions.keys())}}};"
    close_braces_line = "}"
    return '\n'.join([first_line,
                      *known_reals_lines,
                      *ode_metabolite_lines,
                      *derivation_lines,
                      return_line,
                      close_braces_line])


def build_kinetics_function(known_reals,
                            kinetic_functions,
                            ode_metabolites,
                            derived_quantity_expressions):
    first_block = """vector get_kinetics(vector ode_metabolites,
                    vector kinetic_parameters,
                    real[] known_reals){"""

    known_reals_lines = [
        f"  real {kr} = known_reals[{i+1}];"
        for i, kr in enumerate(known_reals.keys())
    ]
    ode_metabolite_lines = [
        f"  real {metabolite} = ode_metabolites[{i+1}];"
        for i, metabolite in enumerate(ode_metabolites)
    ]
    kinetic_parameter_lines = [
        f"  real {parameter} = kinetic_parameters[{i+1}];"
        for i, parameter in enumerate(kinetic_parameters)
    ]
    derived_quantity_calculation_line = (
        f"  real derived_quantities[{len(derived_quantity_expressions.keys())}] = "
        "get_derived_quantities(ode_metabolites, known_reals);"
    )
    derived_quantity_unpack_lines = [
        f"  real {quantity} = derived_quantities[{i+1}];"
        for i, quantity in enumerate(derived_quantity_expressions.keys())
    ]
    kinetic_expression_lines = [
        f"  real {parameter} = {expression}"
        for parameter, expression in kinetic_functions.items()
    ]
    return_line = f"  return [{', '.join(kinetic_functions.keys())}]';"
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


def build_ode_function(odes, kinetic_functions):
    definition_line = "vector get_odes(vector fluxes){"
    flux_unpack_lines = [
        f"  real {flux} = fluxes[{i+1}];"
        for i, flux in enumerate(kinetic_functions.keys())
    ]
    return_line_open = "  return ["
    return_body_lines = [
        f"    {expression}  // {flux}"
        for flux, expression in odes.items()
    ]
    return_body_lines[-1] = return_body_lines[-1].replace(',', '')
    return_close_line = "  ]';"
    close_braces_line = "}"
    return '\n'.join([definition_line,
                      *flux_unpack_lines,
                      return_line_open,
                      *return_body_lines,
                      return_close_line,
                      close_braces_line])


def build_steady_state_function():
    return """vector steady_state_equation(vector ode_metabolites,
                             vector kinetic_parameters,
                             real[] known_reals,
                             int[] known_ints){
    return get_odes(get_kinetics(ode_metabolites, kinetic_parameters, known_reals));\n}"""
            

if __name__ == '__main__':
    # parse input
    metabolites = get_all_metabolites(INPUT_FILE)
    ode_metabolites = get_ode_metabolites(INPUT_FILE)
    kinetic_parameters = get_named_quantity(INPUT_FILE, 'kinetic parameter')
    derived_quantity_expressions = get_derived_quantities(INPUT_FILE)
    kinetic_functions = get_kinetic_functions(INPUT_FILE)
    odes = get_odes(INPUT_FILE)
    known_reals = get_known_reals(INPUT_FILE)

    derived_quantity_function = build_derived_quantity_function(
        known_reals,
        ode_metabolites,
        derived_quantity_expressions
    )
    kinetics_function = build_kinetics_function(
        known_reals,
        kinetic_functions,
        ode_metabolites,
        derived_quantity_expressions,
    )
    ode_function = build_ode_function(odes, kinetic_functions)
    steady_state_function = build_steady_state_function()

    out = '\n\n'.join([derived_quantity_function,
                       kinetics_function,
                       ode_function,
                       steady_state_function])
    output_file = open(OUTPUT_FILE, "w")
    output_file.write(out)
    output_file.close()


