=====================
Specifying input data
=====================

This document explains how to specify kinetic models, prior parameter
distributions and experimental data in a form that Maud can use.

Overview
========

Maud inputs are specified as `toml <https://github.com/toml-lang/toml>`_ files
with three main components: a description of the kinetic model that is being
analysed, a description of some experiments and a specification of prior
distributions that represent pre-experimental information about the network's
kinetic parameters, thermodynamic parameters, unbalanced metabolite
concentrations and enzyme concentrations.

For some working examples of full inputs see `here
<https://github.com/biosustain/Maud/tree/master/data/in>`_.

**NB** The fields depicted in the examples below are all required, except for a
 few optional cases which should be explicitly highlighted.

Specifying a kinetic model
==========================

Kinetic models in Maud input files have three components: compartments,
metabolites and reactions. All of these are specified as tables at the top
level of the input file.

A compartment must have an id, a name and a volume. Here is an example
compartment specification:

.. code:: toml

    [[compartments]]
    id = 'cytosol'
    name = 'cytosol'
    volume = 1

The units for the :code:`volume` field are arbitrary.

A metabolite must have an id, a name, compartment (this should be the id of a
compartment) and a property 'balanced' specifying whether its concentration
should be constant at steady state. Here is an example:

.. code:: toml

    [[metabolites]]
    id = 'AMP'
    name = 'adenosine monophosphate'
    balanced = false
    compartment = 'cytosol'

A reaction can be specified as follows:

.. code:: toml

    [[reactions]]
    id = 'FBA'
    name = 'FBA'
    stoichiometry = { f16p_c = -1, dhap_c = 1, g3p_c = 1 }
    [[reactions.enzymes]]
    id = 'FBA'
    name = 'FBA'
    mechanism = "ordered_unibi"
    allosteric_inhibitors = ['AMP']  # optional

Reaction level information is specified under :code:`[[reactions]]`, and
enzyme-specific information goes under :code:`[[reactions]]`. The stoichiometry
property should map metabolite ids to numerical stoichiometries with arbitrary
units. The mechanism property must be one of the mechanisms that Maud
supports - these can be found in the source code file
`big_k_rate_equations.stan
<https://github.com/biosustain/Maud/blob/master/src/maud/stan_code/big_k_rate_equations.stan>`_. The
optional property allosteric_inhibitors must be a list containing ids of
metabolites that feature in the network.

Specifying experiments
======================

Information about experiments comes in a table called :code:`experiments`,
which can have arbitrarily many entries. Here is an example specification of an
experiment:

.. code:: toml

    [[experiments]]
    id = 'condition_1'
    metadata = "Condition 1"
    metabolite_measurements = [
      { target_id = 'glc__D_c', value = 0.6, uncertainty = 0.1},
      { target_id = 'g6p_c', value = 1.2, uncertainty = 0.1},
      { target_id = 'f6p_c', value = 0.3, uncertainty = 0.1},
      { target_id = 'f16p_c', value = 2.8, uncertainty = 0.1},
      { target_id = 'g3p_c', value = 0.067, uncertainty = 0.1},
      { target_id = 'dhap_c', value = 1.58, uncertainty = 0.1},
      { target_id = '13dpg_c', value = 0.0016, uncertainty = 0.1},
    ]
    reaction_measurements = [
      { target_id = 'GCLt1', value = 1.99, uncertainty = 0.0019},
    ]

Units here are arbitrary, but the values must agree with the rest of the model.

Specifying priors
=================

Priors come in a toml table called :code:`priors`, which must have exactly four
entries: :code:`kinetic_parameters`, :code:`thermodynamic_parameters`
:code:`enzymes` and :code:`unbalanced_metabolites`.

Thermodynamic parameters are specified using this syntax:

.. code:: toml

    [priors.thermodynamic_parameters]
    marginal_dgs = [
      { target_id = 'GLCT1', location = 1, scale = 0.05 },
      { target_id = 'HEX1', location = -17.3, scale = 0.9 },
      { target_id = 'PGI', location = 2.5, scale = 0.8 },
      { target_id = 'PFK', location = -15, scale = 1.3 },
      { target_id = 'FBA', location = 19.8, scale = 1.0 },
      { target_id = 'TPI', location = -5.5, scale = 1.1 },
      { target_id = 'GAPD', location = 7.8, scale = 0.8 },
      { target_id = 'PGK', location = 18.5, scale = 0.9 },
    ]

The :math:`\Delta G` parameters are specified in units of kJ/mol. Each location
and scale input denotes the mean and standard deviation of a normal
distribution over possible values of the :math:`\Delta G` parameter for the
corresponding reaction. These distributions are independent - in future we hope
to implement correlated :math:`\Delta G` priors through separate properties
:code:`mu_dg` and :code:`cov_matrix_dg`.

The :code:`kinetic_parameters` priors should specify marginal kinetic parameter
distributions as follows:

.. code:: toml
    
    [priors.kinetic_parameters]
    GCLt1 = [
      {target_id = 'Kcat1', location = 3.35, scale = 0.1},
      {target_id = 'Ka', location = 0.9, scale = 0.1},
      {target_id = 'Kp', location = 0.9, scale = 0.1},
    ]
    HEX1 = [
      { target_id = 'Kcat1', location = 63.2, scale = 0.1},
      { target_id = 'Ka', location = 0.15, scale = 0.1},
      { target_id = 'Kb', location = 0.293, scale = 0.1},
      { target_id = 'Kp', location = 30, scale = 0.1},
      { target_id = 'Kq', location = 0.23, scale = 0.1},
    ]
    ...

There should be an entry here for every enzyme id in the kinetic model,
containing a line with a :code:`target_id` corresponding to every kinetic
parameter in the enzyme's mechanism.

The kinetic parameters' units are effectively set by those of the :math:`\Delta
G` parameters, through the equality :math:`keq = \exp(\frac{\Delta G}{-RT})`
and the Haldane relationships linking :math:`keq` parameters with other kinetic
parameters.

**NB** Even though kinetic parameters have to be greater than zero and have
lognormal prior distributions, the :code:`location` in these toml inputs are
specified on the standard scale. On the other hand, the :code:`scale` inputs
are interpreted on the log scale with base :math:`e`, representing
multiplicative rather than additive uncertainty.

Priors for steady state enzyme and unbalanced metabolite concentrations are
specified as a series of tables - one for each experiment id - with the
:code:`target_id` inputs corresponding to enzyme ids or metabolite ids. Here is
an example for an input with one experiment called :code:`condition_1`:

.. code:: toml

    [priors.enzymes]
    condition_1 = [
      { target_id = 'GCLt1', location = 1, scale = 0.05 },
      { target_id = 'HEX1', location = 0.062, scale = 0.05 },
      { target_id = 'PGI', location = 0.138, scale = 0.05 },
      { target_id = 'PFK', location = 0.047, scale = 0.05 },
      { target_id = 'FBA', location = 1.34, scale = 0.05 },
      { target_id = 'TPI', location = 0.295, scale = 0.05 },
      { target_id = 'GAPD', location = 0.007, scale = 0.05 },
      { target_id = 'PGK', location = 0.258, scale = 0.05 },
    ]
    
    [priors.unbalanced_metabolites]
    condition_1 = [
      { target_id = 'glc__D_e', location = 10, scale = 1.0 },
      { target_id = 'atp_c', location = 3.95, scale = 0.05 },
      { target_id = 'adp_c', location = 1.72, scale = 0.05 },
      { target_id = 'nad_c', location = 1.41, scale = 0.05 },
      { target_id = 'nadh_c', location = 0.178, scale = 0.05 },
      { target_id = '3pg_c', location = 0.52, scale = 0.05 },
    ]

As with kinetic parameters, the locations are absolute and the scales are
log-scale. The units are arbitrary. When setting them, bear in mind that Stan
tends to work best when most numbers are reasonably close to zero.

