# Script Summary

## export-sbml
Contrary to the name of the script, `export-sbml.py` does not export an sbml, but a yaml
file that is compatible with (yaml2sbml)[https://github.com/yaml2sbml-dev/yaml2sbml]. by using
this script and then the **yaml2sbml** package you will be able to generate an sbml to use in
COPASI.

### Input
* relative_path_Maud_Input
* relative_path_output
* selected_experiment (takes the first one if left `None`)

### Output
* YAML description of system of ODEs with relavant parameter values

### Limitations
* Currently doesn't handle phosphorylation
* Metabolites/Reactions/Drains cannot begin with numbers, even though it's a standard for BiGG identifiers
* Cannot handle the `-` character in any way.
* Doesn't work with priors specified using percentiles