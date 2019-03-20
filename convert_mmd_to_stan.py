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
            

if __name__ == '__main__':
    metabolites = get_all_metabolites(TARGET_FILE)
    ode_metabolites = get_ode_metabolites(TARGET_FILE)
    conserved_totals = get_conserved_totals(TARGET_FILE)
    compartments = get_named_quantity(TARGET_FILE, 'compartment')
    kinetic_parameters = get_named_quantity(TARGET_FILE, 'kinetic parameter')
    global_quantities = get_named_quantity(TARGET_FILE, 'global quantity')
    derived_quantities = get_derived_quantities(TARGET_FILE)
    kinetic_functions = get_kinetic_functions(TARGET_FILE)
    odes = get_odes(TARGET_FILE)

    print(list(metabolites))


