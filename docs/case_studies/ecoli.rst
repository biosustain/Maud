====================================
A simple model of E. coli glycolysis
====================================

This document shows how to use Maud to model the following simple but realistic
biological system:

.. figure:: ecoli_glycolysis.png
    :scale: 50%

This drawing was obtained by deleting the vast majority of reactions from the
`IJO1366 E. coli central metabolism model
<https://escher.github.io/#/app?map=iJO1366.Central%20metabolism&tool=Builder&model=iJO1366>`_
using the online modelling tool `escher <https://escher.github.io/#/>`_.

The target system has 5 reactions, including the interesting
phosphofructokinase reaction, which is thought to be highly regulated and
instrumental for alleviating `redox stress <http://linkinghub.elsevier.com/retrieve/pii/S2405471218301492>`_. The reactions are part of the glycolysis
pathway, which converts glucose into pyruvate and some ATP.

Constructing a suitable input
=============================

In order to model this system with Maud, we first need to decide how to
represent it, and some information about it, in `Maud's input format
<usage/inputting.html>`.

The full input folder can be found in `Maud's GitHub repository
<https://github.com/biosustain/Maud/blob/master/tests/data/ecoli_small>`_. The
following section explains how it was constructed.


Kinetic Model
-------------

The kinetic model is a representation of the reactions involved in the system
defined above using enzyme kin

First, we need to decide how we will define a steady state for this system. In
practical terms this means we need to specify which of the metabolites in the
network we want to treat as "balanced", i.e., such that at steady state their
concentration should be constant. In this case the balanced metabolites we
chose are f6p, fdp and dhap.

The metabolites h2o and h --i.e. water and hydrogen ion--are involved in the
PFK and FBP reactions, but are typically ignored in kinetic analyses. In order
to avoid problems with interpreting prior information from the literature, we
can ignore h by simply leaving it out of our representation of the PFK
reaction. However, we cannot fully ignore the role of h2o in the FBP reaction
because h2o has non-zero formation energy: leaving it out would result in incorrect Gibbs energy estimates. 
We can inform Maud that it needs to make the necessary thermodynamic adjustment by
adding a non-empty `water_stoichiometry` field to the FBP reaciton.

Finally, for the sake of simplicity this case study ignores all regulation,
even though this is not realistic as the PFK reaction is highly regulated.

The kinetic model file therefore looks as follows:

.. code-block:: toml

    metabolites = [
      {id="g6p", name="D-Glucose 6-phosphate", compartment="c", balanced=false},
      {id="f6p", name="D-Fructose 6-phosphate", compartment="c", balanced=true},
      {id="fdp", name="D-Fructose 1,6-bisphosphate", compartment="c", balanced=true},
      {id="adp", name="ADP C10H12N5O10P2", compartment="c", balanced=false},
      {id="atp", name="ATP C10H12N5O13P3", compartment="c", balanced=false},
      {id="pi", name="Phosphate", compartment="c", balanced=false},
      {id="dhap", name="Dihydroxyacetone phosphate", compartment="c", balanced=true},
      {id="g3p", name="Glyceraldehyde 3-phosphate", compartment="c", balanced=true},
    ]

    [[compartments]]
    id = "c"
    name = "cytosol"
    volume = 1

    [[reactions]]
    id = "PGI"
    name = "Glucose-6-phosphate isomerase"
    enzymes = [{id = "PGI", name = "Glucose-6-phosphate isomerase"}]
    stoichiometry = {g6p_c = -1, f6p_c = 1}
    mechanism = "reversible_modular_rate_law"

    [[reactions]]
    id = "PFK"
    name = "Phosphofructokinase"
    enzymes = [{id = "PFK", name = "Phosphofructokinase"}]
    stoichiometry = {atp_c = -1, f6p_c = -1, adp_c = 1, fdp_c = 1}
    mechanism = "irreversible_modular_rate_law"

    [[reactions]]
    id = "FBP"
    name = "Fructose-bisphosphatase"
    water_stoichiometry = 1
    enzymes = [{id = "FBP", name = "Fructose-bisphosphatase"}]
    stoichiometry = {f6p_c = -1, fdp_c = 1, pi_c = -1}
    mechanism = "reversible_modular_rate_law"

    [[reactions]]
    id = "FBA"
    name = "Fructose-bisphosphate aldolase"
    enzymes = [{id = "FBA", name = "Fructose-bisphosphate aldolase"}]
    stoichiometry = {dhap_c = 1, fdp_c = -1, g3p_c = 1}
    mechanism = "reversible_modular_rate_law"

    [[reactions]]
    id = "TPI"
    name = "Triose-phosphate isomerase"
    enzymes = [{id = "TPI", name = "Triose-phosphate isomerase"}]
    stoichiometry = {dhap_c = -1, g3p_c = 1}
    mechanism = "reversible_modular_rate_law"

    [[drains]]
    id = "g3p_drain"
    name = "g3p_drain"
    stoichiometry = { g3p_c = -1 }


Priors
------

Priors for the 8 metabolites' formation energies were found using `equilibrator
<http://equilibrator.weizmann.ac.il/>`_, and are as follows:

.. csv-table:: Formation energy priors

    parameter_type,metabolite_id,mic_id,enzyme_id,drain_id,phos_enz_id,experiment_id,location,scale,pct1,pct99
    dgf,g6p,,,,,,-1336.3,1.3,,
    dgf,f6p,,,,,,-1333.8,1.3,,
    dgf,pi,,,,,,-1073.3,1.5,,
    dgf,adp,,,,,,-1440.8,2.4,,
    dgf,atp,,,,,,-2313.0,3.0,,
    dgf,fdp,,,,,,-2220.9,2.1,,
    dgf,g3p,,,,,,-1106.4,1.3,,
    dgf,dhap,,,,,,-1111.9,1.1,,

This specification highlights a limitation of Maud's prior model: currently
Maud can only specify priors for formation energies as independent normal
distribution. In reality, there is information available not just about the
marginal values of each metabolite's formation energy, but also about
correlations between them. This is because formation energies are typically
estimated based on observations that depend on linear combinations of formation
energies. For example, the formation energies of atp and adp are estimated
using observations of the adenylate kinase reaction; these observations are
determined by a linear combination of the formation energies of atp, adp and
amp. These observations constrain the sum of atp and adp's formation energies
more closely than the marginal values. The result of this limitation is that
Maud's prior model assigns weight to formation energy configurations that are
very unlikely given the underlying information, something that should be fixed
in a future implementation.

Priors for reaction :math:`k_{cat}` and :math:`k_m` parameters are taken from
the `sabio <http://sabio.h-its.org/>`_ database, and are specified in the toml
input as follows:

.. csv-table:: :math:`k_{cat}` and :math:`km` priors

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
    km,,fdp_c,FBP,,,,16.0,0.2,,
    km,,f6p_c,FBP,,,,0.689,1.5,,
    km,,pi_c,FBP,,,,1.0,1.5,,
    km,,fdp_c,FBA,,,,0.02,0.2,,
    km,,g3p_c,FBA,,,,0.03,0.2,,
    km,,dhap_c,FBA,,,,0.13,0.2,,
    km,,dhap_c,TPI,,,,2.16,1.5,,
    km,,g3p_c,TPI,,,,200.0,0.2,,


Experimental data
-----------------

For this case study we pretend that one experiment was carried out, with the
following artificial but approximately realistic results:


.. csv-table:: Experiments

    measurement_type,target_id,experiment_id,measurement,error_scale
    mic,f6p_c,Evo04ptsHIcrrEvo01EP,0.6410029,0.146145
    mic,fdp_c,Evo04ptsHIcrrEvo01EP,4.5428601,0.237197
    mic,dhap_c,Evo04ptsHIcrrEvo01EP,1.895018,0.078636
    mic,f6p_c,Evo04Evo01EP,0.6410029,0.146145
    mic,fdp_c,Evo04Evo01EP,4.5428601,0.237197
    mic,dhap_c,Evo04Evo01EP,1.895018,0.078636
    flux,PGI,Evo04ptsHIcrrEvo01EP,4.08767353555,1
    flux,PGI,Evo04Evo01EP,4.08767353555,1


Fitting the model
=================

We can fit the model from Maud's root directory by running the following
command in a suitable python environment (see `here
<usage/post_installation_usage.html>` for full details about how to run Maud).

.. code-block:: bash

    maud sample tests/data/ecoli_small


Analysing the results
=====================

After a little while, Stan's sampler has finished, an output directory has been
created and populated and Maud has printed the following diagnostic information:

.. code-block:: bash

    Checking sampler transitions treedepth.
    Treedepth satisfactory for all transitions.

    Checking sampler transitions for divergences.
    No divergent transitions found.

    Checking E-BFMI - sampler transitions HMC potential energy.
    E-BFMI satisfactory for all transitions.

    Effective sample size satisfactory.

    Split R-hat values satisfactory all parameters.

The diagnostic message raises no warnings, indicating that Maud's output files
probably represent draws from the posterior distribution defined by our input.

Investigating the marginal posterior distributions for metabolite
concentrations, the results appear broadly plausible.

.. figure:: conc.png

Similarly, the marginal posteriors for reaction fluxes are close to the
measured value of -0.5 for FBP and 4.08 for other reactions:

.. figure:: conc.png

Finally, the marginal posteriors for kinetic parameters are also plausible,
though the :math:`k_{cat}` parameter for the TPI reaction is very high at
around 10000.

.. figure:: kinetic_params.png
