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

Steps to address computational issues
=====================================

This section explains steps that we have taken to make Maud's problem more
computationally tractable.

Use fast equation solvers
-------------------------

Maud currently uses the CVODES `backwards differentiation formula
<http://sundials.wikidot.com/bdf-method>`_ ODE solver via its Stan interface to
calculate steady state metabolite concentrations and fluxes. This solver
`compares favourably with alternatives
<http://www.stochasticlifestyle.com/comparison-differential-equation-solver-suites-matlab-r-julia-python-c-fortran/>`_
in efficiency and is well integrated into Stan.

Better results could possibly be achieved by using a custom solver that is
specifically tailored to metabolic steady state problems. However this would
require not just writing a new solver but also integrating it with Stan's
automatic differentiation framework.

Another approach which we have explored is to calculate steady states using a
hybrid ODE solver / algebra solver approach as outlined in `this paper
<https://zenodo.org/record/1284375>`_. When implemented previously using `Stan's
built in Powell solver
<https://mc-stan.org/docs/2_24/functions-reference/functions-algebraic-solver.html>`_
we achieved mixed results. For simple problems we saw speed increases, but for
more complicated problems the algebra solver would often fail to find a
solution. Since we tried this, Stan has introduced an interface to the Sundials
Newton solver, which should be more robust. We therefore think it is worth
exploring this approach again.


Parallelism
-----------

Maud currently uses two forms of parallelism it runs separate markov chains in
separate processes and calculates likelihood sums using different threads via
Stan's `reduce_sum
<https://mc-stan.org/docs/2_24/stan-users-guide/reduce-sum.html>`_
functionality.

Further options include MPI and GPU parallelism. MPI parallelism can be used
via Stan's `map_rect
<https://mc-stan.org/docs/2_24/stan-users-guide/map-rect.html>`_ functionality
to distribute separable calculations between different processes. For example,
it should be possible to calculate ODE solutions for different experiments in
parallel. This functionality will be implemented when the speed gains can be
shown to justify the cost in complexity.

Finally, `Stan is introducing support for GPU parallelism
<https://arxiv.org/abs/1907.01063>`_. Exploiting these new features should not
require any work on our part beyond configuring cmdstan to use GPU resources
where possible. It may be possible to add new GPU functionality to Stan that
specifically targets Maud's speed bottlenecks, but this option has not yet been
explored.


Efficient Stan coding
---------------------

The following general principles tend to lead to faster Stan programs and have
been kept to where possible:

- performing calculations in the transformed data and generated quantities
  blocks in preference to the transformed parameters and model blocks
- keeping parameters close to unit scale
- using vectorised operations rather than loops
- parameterising the model so that the parameters are as uncorrelated

The final point about less correlated parameterisation is under active
development - it is likely that the thermodynamic component of the model can be
reparameterised so as to reduce parameter correlation.


Better priors
-------------

There are several as yet unexplored ways in which changing Maud represents
prior information could improve performance.

It is possible that the log-normal distribution that Maud uses to express prior
information about non-negative unknowns is not optimal, allocating too much
prior mass to the tails, contrary to both the available information and optimal
computation. It is therefore worth exploring alternatives like the gamma
distribution.

Allowing information about parameter correlations to be expressed might also
improve computaion.

Finally, Maud does not currently support input of prior information about the
values of fluxes or balanced metabolite concentrations, mainly because this is
technically difficult to implement. Both of these features could make Maud's
posterior distributions less degenerate and easier to explore, thereby
improving computation.
