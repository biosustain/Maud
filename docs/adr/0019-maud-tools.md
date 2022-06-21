# 19. Create a separate `maudtools` repository for maud related tools

Date: 2021-03-30

## Status

Accepted

## Context

Maud's core functionality is as follows:

- Specify input format
- Validate input
- Define statistical model in Stan code
- Run the statistical model by connecting the input and Stan code with Stan
- Store Stan output in a nice format

However Maud currently also provides some other functionality:

- Generate inits from a previous run with the command `maud generate-inits`
- Generate a prior template file with the command `maud generate-prior-template`
- Fetch formation energy priors from equilibrator with the script `scripts/get_dgf_priors_from_equilibrator.py`
- Create a `yaml2sbml` compatible yaml file with the script `scripts/export_yaml.py`
- Create plots with the scripts `src/maud/plotting/plot_oos.py` and `src/maud/plotting/plot_samples.py`

This is inconsistent: sometimes you use the main cli to do a non-core task,
sometimes you use a script. The scripts are in different places.

There are more things we would like to provide that are outside the core
functionality:

- More functionality for conveniently creating Maud inputs, e.g. generating a
  kinetic model file from sbml.
- Online monitoring
- More and better plots

## Decision

Create a separate repository for Maud related tools, called `maudtools`. Move
the existing non-core functionality to the new repository and write the new
things we want there.

## Consequences

This will solve Maud's feature creep problem without reducing the number of
things we can do. As a result it will become simpler to navigate and maintain
the Maud repository.

It will also make for more predictable user experience and make it easier to
test and validate the extra functionality.

On the other hand this change will make it necessary to keep the two projects in
sync. This can be made easier by adding proper versioning to Maud and testing
`maudtools` functionality against a specific version.
