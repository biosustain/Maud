===========
Computation
===========

Fitting Bayesian statistical models of large metabolic networks poses
distinctive computational problems. This document explains what the problems
are and how Maud attempts to address them.

Goal
====

Maud aims to fit Bayesian models of metabolic networks with about 100 internal
reactions, each with its kinetics represented by the modular rate law and MWC
formalism as described in the section on :doc:`Maud's statistical model
<../theory/statistics>`, and around 20 experiments.

Inputs
======

There should be informative prior information about all kinetic parameters,
compound formation energies and unbalanced metabolite
concentrations. Compartment volumes, measurements and measurement errors and a
kinetic model definition must be entered exactly.

This input information is put in a `toml <https://github.com/toml-lang/toml>`_
document.

We are currently researching how accurate the measurements and priors should be
in order to produce accurate results.

How Maud draws posterior samples
================================

Maud draws posterior samples using Stan and cmdstanpy. When a user inputs a
:code:`maud sample ...` command, the chosen toml file, and any other
configuration options are read and translated into :doc:`Maud's data model
</data_model>`, then initial parameter values and input data in Stan format are
written as json files. If necessary, the Stan program `<inference_model.stan>`
is then compiled to c++ and then binary. Next, Maud uses cmdstanpy to trigger
Stan's sampler using the specified input data and initial conditions. The
resulting samples are written to csv files in Stan's output format.

Stan's sampler uses an adaptive Hamiltonian Monte Carlo algorithm to
efficiently explore the posterior probability distribution defined by the
:doc:`statistical model <../Theory/Statistics>` And Input Data.


Why Use Stan And Adaptive Hmc?
------------------------------

Our Approach Using Adaptive Hmc Is Preferable To Alternatives Based On
Rejection Sampling Because The Parameter Space That Needs To Be Explored In
Scientifically Interesting Problems Has A Large Number Of Dimensions. Rejection
Sampling Is Known To Be Unreliable In High-Dimensional Problems.

Compared To Other Markov Chain Monte Carlo Approaches, Ours Has The Following
Notable Advantages:

- It Is Highly Efficient
- The well-tested and fast CVODES ODE solver and algebra solver are available
- Software complexity can be kept low compared to a custom sampler
- Effective diagnostic tools are available, so we can be relatively confident,
  in the absence of diagnostic warnings, that Maud's output really represents
  draws from the target posterior distribution.
- It is under active development, allowing Maud to passively benefit from
  upstream work.
  
In future we plan to explore alternative approaches using approximate tools
like variational inference. However, since these can sometimes be unreliable
and lack equivalent diagnostics to adaptive Hamiltonian Monte Carlo

Bottlenecks
===========

Even though it is more efficient than the available alternatives, exploring the
whole parameter space using Stan's adaptive Hamiltonian Monte Carlo algorithm
requires many evaluations of the posterior density and its gradients with
respect to the values of all unknown parameters. This is because, for each
sample, the algorithm needs to calculate forwards and backwards Hamiltonian
trajectories using a numerical leapfrog integration algorithm. These
trajectories must traverse a large distance while avoiding discretisation
errors. With clever adaptation it is possible to transform the posterior
surface so that it can be explored more efficiently, but depending on the
posterior distribution a large variation is possible in the number of
leapfrog steps per iteration.

Maud's problem is unusual in that a large system of equations needs to be
solved in order to evaluate the posterior density and gradients. Maud's
computational cost per leapfrog step is therefore relatively high. The
posterior surfaces that Maud needs Stan to explore are also relatively
complicated due to the many non-linearities in the target systems.

