import re

TARGET_FILE = "Ecoli glycolisis.mmd"


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
    
def get_all_metabolites(target_file):
    with open(target_file, 'r') as f:
        for l in f:
            if 'concentration of metabolite' in l:
                yield re.findall("concentration of metabolite '(.*?)'", l)[0]



def get_ode_metabolites(target_file):
    with open(target_file, 'r') as f:
        for l in f:
            if l[:4] == 'init':
                met = l.split(' ')[1]
                yield met
            

def get_conserved_totals(target_file):
    with open(target_file, 'r') as f:
        for l in f:
            if 'conserved total' in l:
                ct = remove_non_numeric_bits(l.split(' ')[2])
                yield ct


def get_named_quantity(target_file, quantity_type):
    out = dict()
    with open(target_file, 'r') as f:
        for l in f:
            if quantity_type in l:
                split_line = l.split('=')
                name = split_line[0].strip()
                volume = remove_non_numeric_bits(split_line[1])
                out[name] = volume
    return out


def get_derived_quantities(target_file):
    out = dict()
    on = False
    with open(target_file, 'r') as f:
        for l in f:
            if l == ' \n':
                on = False
            if on:
                quantity = l.split(' ')[0]
                expression = re.findall("\= (.*?);", l)[0].replace('\t', '') + ';'
                out[quantity] = expression
            if l == '{Assignment Model Entities: }\n':
                on = True
    return out


def get_kinetic_functions(target_file):
    out = dict()
    with open(target_file, 'r') as f:
        for l in f:
            if 'kinetic function' in l:
                reaction = re.findall("reaction '(.*?)'", l)[0]
                expression = l.split(' ')[2].replace('\t', '')
                out[reaction] = expression
    return out


def get_odes(target_file):
    out = dict()
    with open(target_file, 'r') as f:
        for l in f:
            if l[:4] == 'd/dt':
                reaction = re.findall('d/dt\((.*?)\)', l)[0]
                expression = l.split(' ')[2].replace('\t', '').replace('FunctionFor', '')
                out[reaction] = expression
    return out


def build_derived_quantity_function(conserved_totals,
                                    ode_metabolites,
                                    derived_quantity_expressions):
    first_block = """real[] get_derived_quantities(vector ode_metabolites,
                              real[] global_quantities,
                              real[] conserved_totals){"""
    conserved_total_lines = [
        f"  real ct_{i+1} = conserved_totals[{i+1}];"
        for i, _ in enumerate(conserved_totals)
    ]
    ode_metabolite_lines = [
        f"  real {metabolite} = ode_metabolites[{i+1}];"
        for i, metabolite in enumerate(ode_metabolites)
    ]
    derivation_lines = [
        f"  real {quantity} = {expression}"
        for quantity, expression in derived_quantity_expressions.items()
    ]
    return_line = f"  return {{{', '.join(derived_quantity_expressions.keys())}}};"
    close_braces_line = "}"
    return '\n'.join([first_block,
                      *conserved_total_lines,
                      *ode_metabolite_lines,
                      *derivation_lines,
                      return_line,
                      close_braces_line])

            

if __name__ == '__main__':
    metabolites = get_all_metabolites(TARGET_FILE)
    ode_metabolites = get_ode_metabolites(TARGET_FILE)
    conserved_totals = get_conserved_totals(TARGET_FILE)
    compartments = get_named_quantity(TARGET_FILE, 'compartment')
    kinetic_parameters = get_named_quantity(TARGET_FILE, 'kinetic parameter')
    global_quantities = get_named_quantity(TARGET_FILE, 'global quantity')
    derived_quantity_expressions = get_derived_quantities(TARGET_FILE)
    kinetic_functions = get_kinetic_functions(TARGET_FILE)
    odes = get_odes(TARGET_FILE)

    dqf = build_derived_quantity_function(conserved_totals, ode_metabolites, derived_quantity_expressions)
    print(dqf)


