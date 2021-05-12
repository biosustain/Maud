=====================
Specifying input data
=====================

This document explains how to specify kinetic models, prior parameter
distributions and experimental data in a form that Maud can use.

Overview
========

Maud inputs are structured directories, somewhat inspired by the `PEtab
<https://github.com/PEtab-dev/PEtab>`_ format. A Maud input directory must
contain a `toml <https://github.com/toml-lang/toml>`_ file called
:code:`config.toml` which gives the input a name, configures how Maud will be
run and tells Maud where to find the information it needs.

For some working examples of full inputs see `here
<https://github.com/biosustain/Maud/tree/master/tests/data>`_.


Specifying a configuration file
===============================

The file :code:`config.toml` **must** contain the top-level fields
:code:`name`, :code:`kinetic_model`, :code:`priors` and :code:`experiments`. It
can also optionally include keyword arguments to the method
`cmdstanpy.CmdStanModel.sample
<https://github.com/stan-dev/cmdstanpy/blob/develop/cmdstanpy/model.py>`_ in
the table :code:`cmdstanpy_config` and control parameters for Stan's ODE solver
in the table :code:`ode_config`.

Here is an example configuration file:

.. code:: toml

    name = "linear"
    kinetic_model = "kinetic_model.toml"
    priors = "priors.csv"
    experiments = "experiments.csv"
    likelihood = true

    [cmdstanpy_config]
    iter_warmup = 200
    iter_sampling = 200
    chains = 4
    save_warmup = true

    [ode_config]
    abs_tol = 1e-4
    rel_tol = 1e-4
    max_num_steps = 1e9
    timepoint = 50


This file tells Maud that a file representing a kinetic model can be found at
the relative path :code:`kinetic_model.toml`, and that information about priors
and experiments are at :code:`priors.csv` and :code:`experiments.csv`
respectively.

The line :code:`likelihood = true` tells Maud to take into account the
measurements in :code:`experiments.csv`: in other words, **not** to run in
priors-only mode.

When Maud samples with this input, it will create 4 MCMC chains, each with 200
warmup and 200 sampling iterations, which will all be saved in the output csv
files. the ODE solver will find steady states by simulating for 50 seconds,
with a step limit as well as absolute and relative tolerances.


Specifying a kinetic model
==========================

Kinetic models files are specified in `toml
<https://github.com/toml-lang/toml>`_ files, which have three obligatory top
level tables, namely compartments, metabolites and reactions. In addition,
kinetic models can include tables representing drains and phosphorylation
reactions.

A compartment must have an id, a name and a volume. Here is an example
compartment specification:

.. code:: toml

    [[compartments]]
    id = 'c'
    name = 'cytosol'
    volume = 1

The units for the :code:`volume` field are arbitrary.

A metabolite must have an id, a name, compartment (this should be the id of a
compartment) and a property 'balanced' specifying whether its concentration
should be constant at steady state. Here is an example:

.. code:: toml

    [[metabolites]]
    id = 'amp'
    name = 'adenosine monophosphate'
    balanced = false
    compartment = 'c'

A reaction can be specified as follows:

.. code:: toml

    [[reactions]]
    id = 'FBA'
    name = 'FBA'
    stoichiometry = { f16p_c = -1, dhap_c = 1, g3p_c = 1 }
    [[reactions.enzymes]]
    id = 'FBA'
    name = 'FBA'
    [[reactions.enzymes.modifiers]]
    modifier_type = 'allosteric_activator'
    mic_id = 'amp_c'

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

Files containing information about experimental measurements should be csvs
with the following fields:

- :code:`measurement_type`: one out of these options:
  - :code:`mic`: stands for metabolite-in-compartment, has the form :code:`<metabolite_id>_<compartment_id>`
  - :code:`flux`
  - :code:`enzyme`
- :code:`target_id`: the id of the thing measured
- :code:`experiment_id`: an id corresponding to the experiment
- :code:`measurement`: the measured value
- :code:`error_scale`: a number representing the accuracy of the measurement

Error scales are interpreted as the standard deviation of a normal distribution
for flux measurements, which can be negative, or as scale parameters of
lognormal distributions for concentration and enzyme measurements, as these are
always non-negative.

Here is an example experiment file:

.. code:: csv

    measurement_type,target_id,experiment_id,measurement,error_scale
    mic,f6p_c,Evo04ptsHIcrrEvo01EP,0.6410029,0.146145
    mic,fdp_c,Evo04ptsHIcrrEvo01EP,4.5428601,0.237197
    mic,dhap_c,Evo04ptsHIcrrEvo01EP,1.895018,0.078636
    mic,f6p_c,Evo04Evo01EP,0.6410029,0.146145
    mic,fdp_c,Evo04Evo01EP,4.5428601,0.237197
    mic,dhap_c,Evo04Evo01EP,1.895018,0.078636
    flux,PGI,Evo04ptsHIcrrEvo01EP,4.08767353555,1
    flux,PGI,Evo04Evo01EP,4.08767353555,1

Units here are arbitrary, but the values must agree with the rest of the model.

Specifying priors
=================

Files with information about priors should be csvs with the following fields:

- :code:`parameter_type`: see below for options and corresponding id fields:
- :code:`metabolite_id`
- :code:`mic_id`
- :code:`enzyme_id`
- :code:`drain_id`
- :code:`phos_enz_id`
- :code:`experiment_id`
- :code:`location`
- :code:`scale`
- :code:`pct1`: first percentile of the prior distribution
- :code:`pct99`: 99th percentile of the prior distribution

Each parameter type has specific required id fields, which are as follows:

- :code:`kcat`: :code:`enzyme_id`
- :code:`km`: :code:`enzyme_id` and :code:`mic_id`
- :code:`dgf`: :code:`metabolite_id`
- :code:`ki`: :code:`enzyme_id`
- :code:`conc_enzyme`: :code:`enzyme_id` and :code:`experiment_id`
- :code:`conc_unbalanced`: :code:`mic_id` and :code:`experiment_id`
- :code:`drain`: :code:`drain_id` and :code:`experiment_id`
- :code:`transfer_constant`: :code:`enzyme_id`
- :code:`diss_r`: :code:`enzyme_id` and :code:`mic_id`
- :code:`diss_t`: :code:`enzyme_id` and :code:`mic_id`
- :code:`kcat_phos`: :code:`phos_enz_id`
- :code:`conc_phos`: :code:`phos_enz_id` and :code:`experiment_id`

Information in id fields other than the required ones will be ignored: for
clarity it is best to leave these empty, as in the example below.

Quantitative prior information must be represented either using the
:code:`location` and :code:`scale` fields or else the :code:`pct1` and
:code:`pct99` fields.

Formation energy priors should have units of kJ/mol. The units for kinetic
parameter priors are effectively set by those of the formation energies,
through the equality :math:`keq = \exp(\frac{\Delta G}{-RT})` and the Haldane
relationships linking :math:`keq` parameters with other kinetic parameters.

Below is an example priors file.

.. code:: csv

    parameter_type,metabolite_id,mic_id,enzyme_id,drain_id,phos_enz_id,experiment_id,location,scale,pct1,pct99
    kcat,,,PGI,,,,126.0,0.2,,
    kcat,,,PFK,,,,110.0,0.2,,
    kcat,,,FBP,,,,24.0,0.2,,
    kcat,,,FBA,,,,7.0,0.2,,
    kcat,,,TPI,,,,9000.0,0.2,,
    km,,g6p_c,PGI,,,,3.0,0.2,,
    km,,f6p_c,PGI,,,,0.16,0.2,,
    km,,f6p_c,PFK,,,,0.04,0.2,,
    km,,atp_c,PFK,,,,0.06,0.2,,
    km,,fdp_c,PFK,,,,15,1.5,,
    km,,adp_c,PFK,,,,0.55,1.5,,
    km,,fdp_c,FBP,,,,16.0,0.2,,
    km,,f6p_c,FBP,,,,0.689,1.5,,
    km,,pi_c,FBP,,,,1.0,1.5,,
    km,,fdp_c,FBA,,,,0.02,0.2,,
    km,,g3p_c,FBA,,,,0.03,0.2,,
    km,,dhap_c,FBA,,,,0.13,0.2,,
    km,,dhap_c,TPI,,,,2.16,1.5,,
    km,,g3p_c,TPI,,,,200.0,0.2,,
    dgf,g6p,,,,,,-1336.3,1.3,,
    dgf,f6p,,,,,,-1333.8,1.3,,
    dgf,pi,,,,,,-1073.3,1.5,,
    dgf,adp,,,,,,-1440.8,2.4,,
    dgf,atp,,,,,,-2313.0,3.0,,
    dgf,fdp,,,,,,-2220.9,2.1,,
    dgf,g3p,,,,,,-1106.4,1.3,,
    dgf,dhap,,,,,,-1111.9,1.1,,
    conc_enzyme,,,PGI,,,Evo04ptsHIcrrEvo01EP,0.033875912,0.06,,
    conc_enzyme,,,FBP,,,Evo04ptsHIcrrEvo01EP,0.00592291,0.047,,
    conc_enzyme,,,FBA,,,Evo04ptsHIcrrEvo01EP,0.0702922488972023,0.19,,
    conc_enzyme,,,TPI,,,Evo04ptsHIcrrEvo01EP,0.020866941,0.13,,
    conc_enzyme,,,PFK,,,Evo04ptsHIcrrEvo01EP,0.018055101,0.13,,
    conc_enzyme,,,FBP,,,Evo04Evo01EP,0.00592291,0.047,,
    conc_enzyme,,,FBA,,,Evo04Evo01EP,0.0702922488972023,0.19,,
    conc_enzyme,,,TPI,,,Evo04Evo01EP,0.0198,0.1,,
    conc_enzyme,,,PFK,,,Evo04Evo01EP,0.0185,0.05,,
    conc_unbalanced,,g6p_c,,,,Evo04ptsHIcrrEvo01EP,2.0804108,0.188651,,
    conc_unbalanced,,adp_c,,,,Evo04ptsHIcrrEvo01EP,0.6113649,0.038811,,
    conc_unbalanced,,atp_c,,,,Evo04ptsHIcrrEvo01EP,5.4080032,0.186962,,
    conc_unbalanced,,g6p_c,,,,Evo04Evo01EP,2.0804108,0.188651,,
    conc_unbalanced,,adp_c,,,,Evo04Evo01EP,0.6113649,0.038811,,
    conc_unbalanced,,atp_c,,,,Evo04Evo01EP,5.4080032,0.186962,,
    drain,,,,g3p_drain,,Evo04ptsHIcrrEvo01EP,,,0.3,1.2
    drain,,,,g3p_drain,,Evo04Evo01EP,,,0.3,1.2

