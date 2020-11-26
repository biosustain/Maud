==================
Statistical Model
==================

This document describes Maud from a statistical point of view.

Maud's statistical model separates the information that might be available
about a metabolic network into three different kinds:

- structural information implicit in a kinetic model
- information contained in directly modelled experiments
- information from other sources

The role of the statistical model is to synthesise these different sources of
information, making it possible to say exactly what is known about a metabolic
network after a series of experiments.

More specifically, the statistical model is a joint probability distribution
:math:`\pi: \Theta \times Y\rightarrow [0,1]` that defines the probability
density :math:`\pi(\theta, y)` of any possible combination of unknown
parameters :math:`\theta` and observations :math:`y`. This model is written
explicitly as a Stan program.

Given a kinetic model and a set of observations :math:`y`, Maud uses Stan's
adaptive Hamiltonian Monte Carlo algorithm to draw samples from the posterior
distribution :math:`p(\theta\mid y)`. Each draw contains a configuration of
unknown parameter values. Quantitative questions about the metabolic network
can be answered by interrogating the ensemble of posterior draws.


Structural Information
=====================================================

Maud assumes that the following structural information is known in advance:

- The volume of each compartment and the metabolites it houses
- The network stoichiometry, i.e. the proportions in which each reaction in the
  network creates and destroys metabolites and which enzymes catalyse which
  reactions. This information is encapsulated in a stoichiometric matrix
  :math:`S` with a row for each metabolite-in-compartment and a column for each
  enzyme.
- Which metabolites-in-compartment modify which enzymes and how
- Which metabolites-in-compartments are 'balanced', i.e. their concentrations
  do not change when the system is in a steady state.

We refer to this structural information collectively as a kinetic model. See
[LINK] for a detailed description of these from a scientific point of view.

The kinetic model defines a system of ODEs - one differential equation for each
metabolite-in-compartment - with the rate of change of each
metabolite-in-compartment :math:`m_i` described by the following equation:

.. math::

  \frac{dm_{i}}{dt} = \sum_{j} S_{ij} v_{j}(\theta, \mathbf{m})

When the system is at a steady state, all the balanced
metabolites-in-compartments have unchanging concentrations, so their entries in
the equation above are zero.

For a given set of parameters, enzyme concentrations and initial metabolite
conditions, there should be a single steady state balanced metabolite
concentration and set of fluxes.

The kinetic model's role in Maud's statistical model is to connect latent
parameters - i.e. :math:`\theta` above - with measureable quantities,
i.e. :math:`\mathbf{m}` and :math:`\mathbf{v}`.


Probabilistic Model
===================

Maud aims to implement a Bayesian probabilistic model where the joint distribution 
:math:`\pi(\theta, y)` of unknowns and observations is factored into a measurement 
model or likelihood :math:`\pi(y \mid \theta)` and a prior model :math:`\pi(\theta)`. 
This section explains how each of these components is constructed.

Likelihood
----------
Maud represents information from experiments that measure enzyme concentrations
and metabolite concentrations using the following regression model, where
:math:`y` is the observation and :math:`\hat{y}` is the unobserved true value
of the experimentally measured quantity:

.. math::

   y \sim LogNormal(\log(\hat{y}), \sigma)

Measurements of reaction fluxes use the following similar regression model:

.. math::

   y \sim Normal(\hat{y}, \sigma)
   

In all three cases it is assumed that for each kind of measurement the error
standard deviation :math:`\sigma` is known (though this number is in general
different for each measured quantity in each experiment). It is the user's
responsibility to choose values that reflect the accuracy of the measurement
apparatus.

The log-normal distribution was chosen to represent experimental metabolite and
enzyme concentration measurements because the apparatuses used to measure these
quantities typically have multiplicative errors. In other words, if the
measured value is large, then the associated error is also proportionally
large.

Summary statistics
++++++++++++++++++

It is common for experimental results to be reported in the form of a sample
mean and standard deviation. It is important to note that for non-negative
quantities like metabolite and enzyme concentrations these summary statistics
will generally not be good values use as :math:`y` and :math:`\sigma` above. If
possible, non-summarised measurement results should be used instead.


Relative measurements
+++++++++++++++++++++

In many realistic cases a measurement apparatus will give comparatively
accurate information about the relative concentrations of some metabolites or
enzymes while being comparatively uninformative as to their absolute
values. While Maud currently does not support this kind of measurement, support
is planned and will take the following form.


Priors
------

Information that does not naturally take the form of an experimental
measurement can be expressed in Maud's prior model. Maud allows users to
specify independent log-normal priors for the following quantities:

- enzyme :math:`kcat` parameters
- enzyme/metabolite-in-compartment :math:`km` parameters
- enzyme transfer constants
- enzyme/metabolite-in-compartment inhibition and dissociation constants
- enzyme concentrations
- unbalanced metabolite concentrations

The distinction between balanced and unbalanced metabolites is also found in
the statistical model. Information about the unbalanced metabolites
can be parsed in the form of a prior, however, due to the difficulty of non-linear transformations, 
balanced metabolites are always evaluated as part of the model likelihood.
The distinction between unbalanced and balanced becomes aparent when considering what
the unbalanced metabolites represent, which is a boundary condition. These
define the outcome of systems of differential equations, in the case of Maud
this happens to be balanced metabolite concentrations and fluxes. And, our
knowledge about the state of each condition is only conveyed through priors
on the boundary conditions:
* unbalanced metabolite concentrations,
* enzyme concentrations, 
* kinetic parameters, and,
* drains.


For metabolite formation energies, which can be both negative and positive
numbers, Maud allows users to specific independent normal priors.

Users are encouraged to choose prior locations and scales using the method by
calculating quantiles. Prior information is often easiest to ellicit in the
form of qualitative statements like "it is very unlikely that :math:`kcat_e` is
higher than 6.8 or lower than 0.4". Information in this form naturally
translates into restrictions on the quantiles of the corresponding marginal
prior distribution - for example that the prior mass for the events
:math:`kcat_e > 6.8` and :math:`kcat_e < 0.4` should each be about 1%. The
prior values can then be calculated as roughly :math:`\mu_{kcat_e} = 0.5003`
and :math:`\sigma_{kcat_e} = 0.6089`.

Maud includes convenience functions for working out priors in this way, which
can be used in a python environment as follows:

.. code::

  In [1]: from maud.utils import get_lognormal_parameters_from_quantiles 

  In [2]: get_lognormal_parameters_from_quantiles(0.4, 0.01, 6.8, 0.99)
  Out[2]: (0.5003159401539531, 0.608940170915830)


Information about fluxes and balanced metabolite concentrations
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

It is currently not possible to include non-experimental information about
fluxes and steady-state concentrations of balanced metabolites.

This is due to a technical limitatation. Since fluxes and steady state
metabolite concentrations are calculated from the values of other parameters by
finding the solution to the ODE system, directly setting priors would introduce
a bias without a compensating Jacobian adjustment. We have not found a way to
introduce this Jacobian adjustment, so Maud unfortunately cannot currently
represent this information.


Multivariate priors
+++++++++++++++++++

Sometimes the non-experimental information about two parameters is not
independent. For example, some linear combinations of formation energies are
known within a relatively small range even though the marginal value of each
component of the linear combination is not well known.

In such cases a multivariate distribution is required in order to express the
available information. This functionality is not yet supported, but will be
soon.
