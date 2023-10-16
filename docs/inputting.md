Specifying input data
=====================

This page explains how to create input data for Maud.

## Overview

Maud inputs are structured directories, somewhat inspired by the [PEtab](
https://github.com/PEtab-dev/PEtab) format. A Maud input directory must contain
a [toml](https://github.com/toml-lang/toml) file called `config.toml` which
gives the input a name, configures how Maud will be run and tells Maud where to
find the information it needs. It must also include a file containing a kinetic
model definition, one specifying independent prior distributions, a file with
information about the experimental setup and another one recording the results
of measurements. Finally, the input folder can also optionally include extra
files specifying non-independent priors and a file specifying initial parameter
values for the MCMC sampler.

For some working examples of full inputs see [here](https://github.com/biosustain/Maud/tree/master/maud/data/example_inputs).

## The configuration file


The file `config.toml` **must** contain these top-level fields:

* `name` String naming the input
* `kinetic_model_file` Path to a `toml` file defining a kinetic model
* `priors_file` Path to a `toml` file specifying priors
* `experiments_file` Path to a `toml` file with information about experiments
* `likelihood` Boolean representing whether to use information from measurements

The following optional fields can also be specified:

* `reject_non_steady` Boolean saying whether to reject draws that enter non-steady states
* `penalize_non_steady` Boolean saying whether to penalize steady state deviations in the likelihood. It cannot be `True` when `reject_non_steady` is `True`.
* `ode_solver_config` Table of configuration options for Stan's ode solver
* `algebra_solver_config` Table of configuration options for Stan's algebra solver
* `cmdstanpy_config` Table of keyword arguments to the cmdstanpy method [`CmdStanModel.sample](https://cmdstanpy.readthedocs.io/en/v1.1.0/api.html#cmdstanpy.CmdStanModel.sample)
* `cmdstanpy_config_predict` Table of overriding sample keyword argments for predictions
* `stanc_options` Table of valid choices for [CmdStanModel](https://cmdstanpy.readthedocs.io/en/v1.1.0/api.html#cmdstanpy.CmdStanModel) argument `stanc_options`
* `cpp_options` Table of valid choices for  [`CmdStanModel](https://cmdstanpy.readthedocs.io/en/v1.1.0/api.html#cmdstanpy.CmdStanModel) argument `cpp_options`
* `variational_options` Arguments for the cmdstanpy method [`CmdStanModel.variational](https://cmdstanpy.readthedocs.io/en/v1.1.0/api.html#cmdstanpy.CmdStanModel.variational)
* `user_inits_file` path to a toml file of initial values
* `steady_state_threshold_abs` absolute threshold for Sv=0 be at steady state
* `steady_state_threshold_rel` relative threshold for Sv=0 be at steady state
* `steady_state_penalty_sigma` standard deviation to penalize deviations from steady state in the likelihood
* `default_initial_concentration` default initial concentration for unmeasured species
* `drain_small_conc_corrector` number for correcting small conc drains
* `molecule_unit` A unit for counting molecules, like 'mol' or 'mmol'
* `volume_unit` A unit for measuring volume, like 'L'
* `energy_unit` A unit for measuring energy, like 'J' or 'kJ'


Here is an example configuration file:

.. code:: toml

    name = "linear"
    kinetic_model_file = "kinetic_model.toml"
    priors_file = "priors.toml"
    experiments_file = "experiments.toml"
    likelihood = true
    steady_state_threshold_abs = 1e-6

    [cmdstanpy_config]
    refresh = 1
    iter_warmup = 200
    iter_sampling = 200
    chains = 4
    save_warmup = true

    [ode_solver_config]
    abs_tol = 1e-4
    rel_tol = 1e-4
    max_num_steps = 1e6

    [algebra_solver_config]
    abs_tol = 1e-4
    rel_tol = 1e-4
    max_num_steps = 1e6

This file tells Maud that a file representing a kinetic model can be found at
the relative path `kinetic_model.toml`, and that priors and experimental
information can be found at `priors.toml`, `experiments.toml`
respectively.

The line `likelihood = true` tells Maud to take into account the
measurements in `experiments.toml`: in other words, **not** to run in
priors-only mode.

When Maud samples with this input, it will create 4 MCMC chains, each with 200
warmup and 200 sampling iterations, which will all be saved in the output csv
files. the ODE solver will find steady states by simulating for 1000 seconds,
with a step limit as well as absolute and relative tolerances.

## The kinetic model file

A Maud input should use exactly one kinetic model file, which is written in the
[toml](https://github.com/toml-lang/toml) markup language and pointed to by the
`kinetic_model` field of the input's `config.toml` file. This section explains
how to write this kind of file.

If it doesn't make sense, make sure to check the
[code](https://github.com/biosustain/Maud/blob/master/maud/data_model/kinetic_model.py)
that tells Maud what a kinetic model should look like.

### name

This top level field is a string describing the kinetic model.

### compartment

A table with the following obligatory fields:

* `id` A string identifying the compartment without any underscore characters.
* `name` A string describing the compartment
* `volume` A float specifying the compartment's volume

Here is an example compartment table:


```toml
compartment = [
  {id = 'c', name = 'cytosol', volume = 1},
  {id = 'e', name = 'external', volume = 1},
]
```

### metabolite

A table with the following obligatory fields:

* `id` A string identifying the metabolite without any underscore characters.
* `name` A string describing the metabolite

Here is an example metabolite table:

```toml
metabolite = [
  {id = "M1", name = "Metabolite number 1"},
  {id = "M2", name = "Metabolite number 2"},
]
```

### metabolite_in_compartment

A table that specifies which metabolites exist in which compartments, and
whether they should be considered balanced or not. The fields in this table are
as follows:

* `metabolite_id` The id of an entry in the `metabolite` table
* `compartment_id` The id of an entry in the `compartment` table
* `balanced` A boolean

For a `metabolite_in_compartment` to be balanced means that its
concentration does not change when the system is in a steady state. Often
metabolites in the external compartment will be unbalanced.

Here is an example `metabolite_in_compartment` table:


```toml
metabolite_in_compartment = [
  {metabolite_id = "M1", compartment_id = "e", balanced = false},
  {metabolite_id = "M1", compartment_id = "c", balanced = true},
  {metabolite_id = "M2", compartment_id = "c", balanced = true},
  {metabolite_id = "M2", compartment_id = "e", balanced = false},
]
```

### enzyme

A table with the following obligatory fields:

* `id` A string identifying the enzyme without any underscore characters.
* `name` A string describing the enzyme
* `subunits` An integer specifying how many subunits the enzyme has.


```toml
    enzyme = [
      {id = "r1", name = "r1ase", subunits = 1},
      {id = "r2", name = "r2ase", subunits = 1},
      {id = "r3", name = "r3ase", subunits = 1},
    ]
```

### reaction

A table with the following obligatory fields:

* `id` A string identifying the reaction without any underscore characters.
* `name` A string describing the reaction
* `mechanism` A string specifying the reaction's mechanism
* `stoichiometry` A mapping representing the stoichiometric coefficient
  for each `metabolite_in_compartment` that the reaction creates or
  destroys.

In addition the following optional fields can be specified:

* `water_stoichiometry` A float indicating the reaction's water stoichiometry
* `transported_charge` A float indicating the reaction's transported charge. If not specified a default value of zero is used, indicating a non-transport reaction or a transport reaction that does not transport any charge.

Valid options for the `mechanism` field are:

* `reversible_michaelis_menten`
* `irreversible_michaelis_menten`
* `drain`

Each key in the `stoichiometry` should identify an existing
`metabolite_in_compartment` using a `metabolite` id and a
`compartment` id, separated by an underscore.

Here is an example of an entry in a reaction table:


```toml
    [[reaction]]
    id = "r1"
    name = "Reaction number 1"
    mechanism = "reversible_michaelis_menten"
    stoichiometry = { M1_e = -1, M1_c = 1}
```

### enzyme_reaction

A table indicating which enzymes catalyse which reactions, with the following fields:

* `enzyme_id` The id of an entry in the `enzyme` table
* `reaction_id` The id of an entry in the `reaction` table

Here is an example `enzyme_reaction` table:

```toml
    enzyme_reaction = [
      {enzyme_id = "r1", reaction_id = "r1"},
      {enzyme_id = "r2", reaction_id = "r2"},
      {enzyme_id = "r3", reaction_id = "r3"},
    ]
```

### allostery

An optional table with the following fields:

* `enzyme_id` The id of an entry in the `enzyme` table
* `metabolite_id` The id of an entry in the `metabolite` table
* `compartment_id` The id of an entry in the `compartment` table
* `modification_type` A string specifying the kind of modification

Valid options for the `modification_type` field are:

* `activation`
* `inhibition`

Here is an example of an entry in a allostery table:


```toml
    [[allostery]]
    enzyme_id = "r1"
    metabolite_id = "M2"
    compartment_id = "c"
    modification_type = "activation"
```

### competitive_inhibition

An optional table with the following fields:

* `enzyme_id` The id of an entry in the `enzyme` table
* `reaction_id` The id of an entry in the `reaction` table
* `metabolite_id` The id of an entry in the `metabolite` table
* `compartment_id` The id of an entry in the `compartment` table


Here is an example of an entry in a allostery table:


```toml
[[competitive_inhibition]]
enzyme_id = "r2"
reaction_id = "r2"
metabolite_id = "M1"
compartment_id = "c"
```

### phosphorylation

An optional table with the following fields:

* `enzyme_id` The id of an entry in the `enzyme` table
* `modification_type` A string specifying the kind of modification

Valid options for the `modification_type` field are:

* `activation`
* `inhibition`

Here is an example of an entry in a allostery table:


```toml
[[phosphorylation]]
enzyme_id = "r1"
modification_type = "activation"
```

## The experiments file

This is a file written in toml, giving information about the input's
experiments, including qualitative information like whether an enzyme was
knocked out, as well as quantitative information like temperature or the results
of measurements.

This section describes this file's fields.

### experiment

An obligatory table containing information that is specific to each of the
input's experiments. All information in an experiments file should belong to an
experiment.

* `id` A string identifying the experiment, without any underscores
* `is_train` A boolean indicating whether to include the experiment in the
  training dataset
* `is_test` A boolean indicating whether to include the experiment in the
  test dataset
* `temperature` A float specifying the experiment's temperature.
* `enzyme_knockouts`: An optional table describing enzyme knockouts. Each
  entry has one field called "enzyme".
* `pme_knockouts`: An optional table describing knockouts of
  phosphorylation modifying enzymes. Each entry has one field called "pme".
* `measurements`: An optional table describing measurements
* `initial_state`: An optional table describing the initial concentrations of balanced metabolites

The measurement table has these fields:

* `target_type` A string specifying what kind of thing was measured:
  either "mic", "flux" or "enzyme".
* `metabolite` A string identifying the metabolite that was measured,
  required if the target type is "mic".
* `compartment` A string identifying the compartment that was measured,
  required if the target type is "mic".
* `enzyme` A string identifying the enzyme that was measured, required if
  the target type is "enzyme".
* `reaction` A string identifying the reaction that was measured, required
  if the target type is "flux".
- `value` The measured value, as a float.
- `error_scale` The measurement error, as a float.

`error_scale` is the standard deviation of a normal distribution for flux
measurements or the scale parameter of a lognormal distribution for
concentration measurements.

The `initial_state` table has these fields:

* `metabolite` A string identifying the metabolite.
* `compartment` A string identifying the compartment.
* `value` The initial state, as a float.

A default initial concentration for all metabolites and experiments can be
specified in the general configuration. In addition, this `initial_state`
table can be specified for a more granular control of the initial vector of
concentrations to solve the steady state ODE. This may be useful when the
concentrations in the model span over several orders of magnitudes and thus a
unique initial value is suboptimal.

## The priors file

This is a toml file for representing non-experimental quantitative information.

You can specify the following parameters:

* `dgf` (can be negative, identified by metabolite)
* `km` (non-negative, identified by metabolite, compartment, enzyme and reaction)
* `kcat` (non-negative, identified by enzyme and reaction)
* `kcat_pme` (non-negative, identified by phosphorylation modifying enzyme)
* `ki` (non-negative, identified by metabolite, compartment, enzyme and reaction)
* `dissociation_constant` (non-negative, identified by metabolite, compartment and enzyme)
* `transfer_constant` (non-negative, identified by enzyme)
* `psi` (non-negative, identified by metabolite, compartment, enzyme and reaction)
* `conc_unbalanced` (non-negative, identified by metabolite, compartment and experiment)
* `drain` (can be negative, identified by reaction and experiment)
* `conc_enzyme` (non-negative, identified by enzyme and experiment)
* `conc_pme` (non-negative, identified by phosphorylation modifying enzyme and experiment)

All priors are optional. If you do not specify a prior for a parameter, then
Maud will use a default prior instead.

To specify an independent prior distribution for a parameter, create a top level table for it
in your priors file, with each entry containing enough information to identify
the target, as well as its prior distribution.

Independent prior distributions are log-normal for non-negative parameters or
normal for possibly-negative ones. A distribution can be identified either by
specifying its location and scale parameters (fields `location` and
`scale`) or its 1% and 99% quantiles (fields `pct1` and
`pct99`). In addition, for non-negative parameters it is possible to use
the field `exploc` to specify the expontial of the location parameter
instead of directly setting the location, which is sometimes easier to
interpret.

### Multivariate priors for formation energy parameters

For most model parameters, it is safe to model the pre-experimental information
as independent. For example, knowing the value of one enzyme's :math:`kcat`
parameter does not significantly narrow down another enzyme's :math:`kcat`
parameter. Thus in this case, and most others, specifying each parameter's
marginal prior distribution is practically equivalent to specifying the full
joint distribution.

However, the available information about formation energy parameters is
typically not independent. In this case the available information is mostly
derived from measurements of the equilibrium constants of chemical
reactions. Knowing the formation energy of one metabolite is often highly
informative as to the formation energy of another metabolite which produced or
destroyed by the same measured chemical reaction. Metabolites with common
chemical groups are also likely to have similar formation energies, introducing
further non-independence.

In some cases this dependence is not practically important, and Maud will work
well enough with independent priors as described above. For other cases, Maud
allows non-independent prior information to be specified in the form of the mean
vector and covariance matrix of a multivariate normal distribution. This
information is specified using a different format, with fields `ids`,
`mean_vector` and `covariance_matrix`, as below:

```toml
dgf = {
  ids = ["M1", "M2"],
  mean_vector = [-1, 2],
  covariance_matrix = [[1, 0], [0, 1]],
}
```

## The initial parameter values file

Initial parameter values can be entered in a `toml` file. This file should
have a table for each parameter whose inits you would like to set, with
identifiers specified in the same way as priors, and initial values specified
using the subfield `init`. For example:

```toml
kcat = [
  {enzyme = "AHC1", reaction = "AHC", init = 234.284},
  {enzyme = "BHMT1", reaction = "BHMT", init = 13.7676},
  {enzyme = "CBS1", reaction = "CBS", init = 7.02307},
  {enzyme = "GNMT1", reaction = "GNMT", init = 10.5307},
  {enzyme = "MAT1", reaction = "METAT", init = 7.89577},
  {enzyme = "MAT3", reaction = "METAT", init = 19.9215},
  {enzyme = "METH-Gen", reaction = "METH", init = 1.15777},
  {enzyme = "MS1", reaction = "MS", init = 1.77471},
  {enzyme = "MTHFR1", reaction = "MTHFR", init = 3.1654},
  {enzyme = "PROT1", reaction = "PROT", init = 0.264744},
]
```

## Fixing parameter values

Sometimes it is useful to model a parameter as if its value were known exactly.
For example, the formation energies of metabolites that are only involved in
irreversible reactions are typically not identified by Maud's statistical
model. Fixing these unidentified parameters can improve computational
performance without loss of any information. Another use might be to test
whether a particular parameter value is feasible at all.

Maud currently supports fixing the values of `dgf` parameters. This can be done
for independent parameter inputs by adding a `fixed_value` field to the
relevant items as shown here:

```toml
dgf = [
  {location = -1, metabolite = "M1", scale = 0.05, fixed_value = -1},
]
```

If a multivariate input is used for `dgf`, fixed values can be specified using
a separate table, like this:

```toml
[dgf]
ids = ["M1", "M2", "M3"]
mean_vector = [1, 2, 3]
covariance_matrix = [[1, 2, 3], [1, 2, 1], [3, 2, 1]]
fixed_values = {"M1" = -1}
```
