=====================
Specifying input data
=====================

This document explains how to create input data for Maud.

.. contents::
   :depth: 2

Overview
========

Maud inputs are structured directories, somewhat inspired by the `PEtab
<https://github.com/PEtab-dev/PEtab>`_ format. A Maud input directory must
contain a `toml <https://github.com/toml-lang/toml>`_ file called
:code:`config.toml` which gives the input a name, configures how Maud will be
run and tells Maud where to find the information it needs. It must also include
a file containing a kinetic model definition, one specifying independent prior
distributions, a file with information about the experimental setup and another
one recording the results of measurements. Finally, the input folder can also
optionally include extra files specifying non-independent priors and a file
specifying initial parameter values for the MCMC sampler.

For some working examples of full inputs see `here
<https://github.com/biosustain/Maud/tree/master/tests/data>`_.


The configuration file
======================

The file :code:`config.toml` **must** contain these top-level fields:

- :code:`name` String naming the input
- :code:`kinetic_model_file` Path to a :code:`toml` file defining a kinetic model
- :code:`priors_file` Path to a :code:`csv` file specifying independent priors
- :code:`experimental_setup_file` Path to a :code:`toml` file specifying the experimental setup
- :code:`measurements_file` Path to a code:`csv` file specifying measurements
- :code:`likelihood` Boolean representing whether to use information from measurements

The following optional fields can also be specified:

- :code:`reject_non_steady` Boolean saying whether to reject draws that enter non-steady states
- :code:`ode_config` Table of configuration options for Stan's ode solver
- :code:`cmdstanpy_config` Table of keyword arguments to the cmdstanpy method `sample <https://cmdstanpy.readthedocs.io/en/v1.0.1/api.html#cmdstanpy.CmdStanModel.sample>`_
- :code:`cmdstanpy_config_predict` Table of overriding sample keyword argments for predictions
- :code:`stanc_options` Table of valid choices for `CmdStanModel <https://cmdstanpy.readthedocs.io/en/v1.0.1/api.html#cmdstanpy.CmdStanModel>`_ argument `stanc_options`
- :code:`cpp_options` Table of valid choices for  `CmdStanModel <https://cmdstanpy.readthedocs.io/en/v1.0.1/api.html#cmdstanpy.CmdStanModel>`_ argument `cpp_options`
- :code:`variational_options` Arguments for CmdStanModel.variational
- :code:`user_inits_file` path to a csv file of initial values
- :code:`dgf_mean_file` path to a csv file of formation energy means
- :code:`dgf_covariance_file` path to a csv file of formation energy covariances
- :code:`steady_state_threshold_abs` absolute threshold for Sv=0 be at steady state
- :code:`steady_state_threshold_rel` relative threshold for Sv=0 be at steady state
- :code:`drain_small_conc_corrector` number for correcting small conc drains

Here is an example configuration file:

.. code:: toml

    name = "linear"
    kinetic_model_file = "kinetic_model.toml"
    priors_file = "priors.csv"
    measurements_file = "measurements.csv"
    experimental_setup_file = "experimental_setup.toml"
    likelihood = true
    steady_state_threshold_abs = 1e-6

    [cmdstanpy_config]
    refresh = 1
    iter_warmup = 200
    iter_sampling = 200
    chains = 4
    save_warmup = true

    [ode_config]
    abs_tol = 1e-4
    rel_tol = 1e-4
    max_num_steps = 1e6
    timepoint = 1e3

This file tells Maud that a file representing a kinetic model can be found at
the relative path :code:`kinetic_model.toml`, and that priors, experimental
setup information and measurements can be found at :code:`priors.csv`,
:code:`experimental_setup.toml` and :code:`measurements.csv` respectively.

The line :code:`likelihood = true` tells Maud to take into account the
measurements in :code:`measurements.csv`: in other words, **not** to run in
priors-only mode.

When Maud samples with this input, it will create 4 MCMC chains, each with 200
warmup and 200 sampling iterations, which will all be saved in the output csv
files. the ODE solver will find steady states by simulating for 1000 seconds,
with a step limit as well as absolute and relative tolerances.

The kinetic model file
======================

A Maud input should use exactly one kinetic model file, which is written in the
`toml <https://github.com/toml-lang/toml>`_ markup language and pointed to by
the :code:`kinetic_model` field of the input's :code:`config.toml` file. This
section explains how to write this kind of file.

If it doesn't make sense, make sure to check the `code that tells Maud what a
kinetic model should look like
<https://github.com/biosustain/Maud/blob/master/src/maud/data_model/kinetic_model.py>`_.

name
----
This top level field is a string describing the kinetic model.

compartment
-----------
A table with the following obligatory fields:

- :code:`id` A string identifying the compartment without any underscore characters.
- :code:`name` A string describing the compartment
- :code:`volume` A float specifying the compartment's volume

Here is an example compartment table:

.. code:: toml

    compartment = [
      {id = 'c', name = 'cytosol', volume = 1},
      {id = 'e', name = 'external', volume = 1},
    ]

metabolite
----------
A table with the following obligatory fields:

- :code:`id` A string identifying the metabolite without any underscore characters.
- :code:`name` A string describing the metabolite

Here is an example metabolite table:

.. code:: toml

    metabolite = [
      {id = "M1", name = "Metabolite number 1"},
      {id = "M2", name = "Metabolite number 2"},
    ]

metabolite_in_compartment
-------------------------

A table that specifies which metabolites exist in which compartments, and
whether they should be considered balanced or not. The fields in this table are
as follows:

- :code:`metabolite_id` The id of an entry in the :code:`metabolite` table
- :code:`compartment_id` The id of an entry in the :code:`compartment` table
- :code:`balanced` A boolean

For a :code:`metabolite_in_compartment` to be balanced means that its
concentration does not change when the system is in a steady state. Often
metabolites in the external compartment will be unbalanced.
  
Here is an example :code:`metabolite_in_compartment` table:

.. code:: toml

    metabolite_in_compartment = [
      {metabolite_id = "M1", compartment_id = "e", balanced = false},
      {metabolite_id = "M1", compartment_id = "c", balanced = true},
      {metabolite_id = "M2", compartment_id = "c", balanced = true},
      {metabolite_id = "M2", compartment_id = "e", balanced = false},
    ]

enzyme
------

A table with the following obligatory fields:

- :code:`id` A string identifying the enzyme without any underscore characters.
- :code:`name` A string describing the enzyme
- :code:`subunits` An integer specifying how many subunits the enzyme has.

.. code:: toml

    enzyme = [
      {id = "r1", name = "r1ase", subunits = 1},
      {id = "r2", name = "r2ase", subunits = 1},
      {id = "r3", name = "r3ase", subunits = 1},
    ]

reaction
--------

A table with the following obligatory fields:

- :code:`id` A string identifying the reaction without any underscore characters.
- :code:`name` A string describing the reaction
- :code:`mechanism` A string specifying the reaction's mechanism
- :code:`stoichiometry` A mapping representing the stoichiometric coefficient
  for each :code:`metabolite_in_compartment` that the reaction creates or
  destroys.

In addition the following optional fields can be specified:

- :code:`water_stoichiometry` A float indicating the reaction's water stoichiometry
- :code:`transported_charge` A float indicating the reaction's transported charge

Valid options for the :code:`mechanism` field are:

- :code:`reversible_michaelis_menten`
- :code:`irreversible_michaelis_menten`
- :code:`drain`

Each key in the :code:`stoichiometry` should identify an existing
:code:`metabolite_in_compartment` using a :code:`metabolite` id and a
:code:`compartment` id, separated by an underscore.

Here is an example of an entry in a reaction table:

.. code:: toml

    [[reaction]]
    id = "r1"
    name = "Reaction number 1"
    mechanism = "reversible_michaelis_menten"
    stoichiometry = { M1_e = -1, M1_c = 1}

enzyme_reaction
---------------

A table indicating which enzymes catalyse which reactions, with the following fields:

- :code:`enzyme_id` The id of an entry in the :code:`enzyme` table
- :code:`reaction_id` The id of an entry in the :code:`reaction` table

Here is an example :code:`enzyme_reaction` table:

enzyme_reaction = [
  {enzyme_id = "r1", reaction_id = "r1"},
  {enzyme_id = "r2", reaction_id = "r2"},
  {enzyme_id = "r3", reaction_id = "r3"},
]

allostery
---------

An optional table with the following fields:

- :code:`enzyme_id` The id of an entry in the :code:`enzyme` table
- :code:`metabolite_id` The id of an entry in the :code:`metabolite` table
- :code:`compartment_id` The id of an entry in the :code:`compartment` table
- :code:`modification_type` A string specifying the kind of modification

Valid options for the :code:`modification_type` field are:

- :code:`activation`
- :code:`inhibition`

Here is an example of an entry in a allostery table:

.. code:: toml

    [[allostery]]
    enzyme_id = "r1"
    metabolite_id = "M2"
    compartment_id = "c"
    modification_type = "activation"

competitive_inhibition
----------------------

An optional table with the following fields:

- :code:`enzyme_id` The id of an entry in the :code:`enzyme` table
- :code:`reaction_id` The id of an entry in the :code:`reaction` table
- :code:`metabolite_id` The id of an entry in the :code:`metabolite` table
- :code:`compartment_id` The id of an entry in the :code:`compartment` table


Here is an example of an entry in a allostery table:

.. code:: toml

    [[competitive_inhibition]]
    enzyme_id = "r2"
    reaction_id = "r2"
    metabolite_id = "M1"
    compartment_id = "c"

phosphorylation
---------------

An optional table with the following fields:

- :code:`enzyme_id` The id of an entry in the :code:`enzyme` table
- :code:`modification_type` A string specifying the kind of modification

Valid options for the :code:`modification_type` field are:

- :code:`activation`
- :code:`inhibition`

Here is an example of an entry in a allostery table:

.. code:: toml

    [[phosphorylation]]
    enzyme_id = "r1"
    modification_type = "activation"

The experimental setup file
===========================

This is a file written in toml, giving qualititative information about the
input's experimental setup.

This section describes this file's fields.

experiment
----------

An obligatory table containing information that is specific to each of the
input's experiments, with the following fields:

- :code:`id` A string identifying the experiment, without any underscores
- :code:`is_train` A boolean indicating whether to include the experiment in the
  training dataset
- :code:`is_test` A boolean indicating whether to include the experiment in the
  test dataset
- :code:`temperature` A float specifying the experiment's temperature.

enzyme_knockout
---------------

An optional table specifying knockouts of enzymes, with the following fields:

- :code:`experiment_id` Id of the knockout's experiment
- :code:`enzyme_id` Id of the enzyme that was knocked out

phosphorylation_knockout
------------------------

An optional table specifying knockouts of phosphorylation effects, with the
following fields:

- :code:`experiment_id` Id of the knockout's experiment
- :code:`enzyme_id` Id of the enzyme whose phosphorylation was knocked out


The measurements file
=====================

This is a csv file with the following fields:

- :code:`measurement_type` A string specifying what kind of thing was measured 
- :code:`target_id` A string identifying the thing that was measured 
- :code:`experiment` A string specifying the measurement's experiment 
- :code:`measurement` The measured value, as a float
- :code:`error_scale` The measurement error, as a float

Valid options for the :code:`measurement_type` field are:

- :code:`mic` Concentration of a :code:`metabolite_in_compartment`
- :code:`enzyme` Concentration of an enzyme
- :code:`flux` Flux of a reaction

`error_scale` is the standard deviation of a normal distribution for flux
measurements or the scale parameter of a lognormal distribution for
concentration measurements.


The priors file
===============

This is a csv file representing pre-experimental information that can be
represented by independent probability distributions. All of the parameters in
Maud's statistical model have either one or two dimensions: priors for one
dimensional parameters are identified by a :code:`row_id`, while two dimensional
parameters require an additional :code:`column_id`.

The priors table has the following fields:

- :code:`parameter` String identifying a parameter
- :code:`row_id` String identifier
- :code:`column_id` String identifier
- :code:`location` Float specifying a location
- :code:`scale` Float specifying a scale
- :code:`pct1`: First percentile of the prior distribution
- :code:`pct99`: 99th percentile of the prior distribution

See the :code:`id_components` fields in the `corresponding code file
<https://github.com/biosustain/Maud/tree/master/src/maud/data_model/stan_variable_set.py>`_
for the correct :code:`row` for each kind of prior. :code:`column_ids` for two
dimensional parameters are always experiment ids.

Prior distributions can either be specified by a location and scale or by a 1st
and 99th percentile, but not both.

Multivariate priors for formation energy parameters
===================================================

The use of a single csv file for priors was motivated by the fact that, for
most model parameters, it is safe to model the pre-experimental information as
independent. For example, knowing the value of one enzyme's :math:`kcat`
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
well enough with independent priors in a csv file as above. For other cases,
Maud allows non-independent prior information to be specified in the form of
the mean vector and covariance matrix of a multivariate normal
distribution. This information is specified as follows.

First, to indicate where to find the required vector and matrix, the fields
:code:`dgf_mean_file` and :code:`dgf_covariance_file` should be added to the
top level of the file :code:`config.toml` in the input folder. For example:

.. code:: toml

    name = "methionine_cycle"
    kinetic_model = "methionine_cycle.toml"
    priors = "priors.csv"
    experiments = "experiments.csv"
    dgf_mean_file = "dgf_prior_mean.csv"
    dgf_covariance_file = "dgf_prior_covariance.csv"

These fields should be paths from the root of the input folder to csv
files. The :code:`dgf_mean_file` should have columns caled :code:`metabolite`
and :code:`prior_mean_dgf`, with the former consisting of ids that agree with
the rest of the input folder (in particular the kinetic model file) and the
latter of non-null real numbers. For example

.. csv-table::

    metabolite,prior_mean_dgf
    5mthf,-778.2999561
    adn,-190.9913035
    ahcys,-330.3885785
    amet,-347.1029509
    atp,-2811.578332
    cyst-L,-656.8334114
    ...

The :code:`dgf_covariance_file` should be a valid covariance matrix surrounded
by metabolite ids. The first column should be called :code:`metabolite` and
populated with ids that are consistent with the other inputs. Subsequent
columns should have names that match the first column. Here is (the start of)
an example:

.. csv-table::

    metabolite,5mthf,adn,ahcys,amet,atp,cyst-L
    5mthf,457895.226,0.023993053,2.911539829,38.09225442,0.023892737,0.610913519
    adn,0.023993053,2.081489779,1.034504533,1.00E-10,0.444288943,0
    ahcys,2.911539829,1.034504533,16.2459485,4.297104388,0.341195482,13.08072127
    amet,38.09225442,1.00E-10,4.297104388,1000025.576,-1.00E-10,2.066261457
    atp,0.023892737,0.444288943,0.341195482,-1.00E-10,2.22005692,0
    cyst-L,0.610913519,0,13.08072127,2.066261457,0,16.61784088
    ...


The initial parameter values file
=================================

Initial parameter values can be entered in a :code:`json` file. This file should
be a valid option for the :code:`inits` argument of the cmdstan sample method.


