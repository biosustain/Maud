==============
Thermodynamics
==============

Each enzyme in a reaction network model has a thermodynamic parameter called
:math:`k_{eq}`, which represents how much energy the corresponding chemical
reaction stores (:math:`k_{eq}` less than 1) or releases (:math:`k_{eq}`
greater than 1) from its environment. This document explains how the laws of
thermodynamics constrain these parameters and how Maud ensures that these
constraints are satisfied.

Why the laws of thermodynamics impose constraints
=================================================

Thermodynamic parameters are constrained in two main ways.

First, each thermodynamic parameter must agree with the kinetic parameters
describing its enzyme, according to the Haldane relationships governing the
enzyme's mechanism. For example, the Haldane relationships for an enzyme with
an ordered unibi mechanism are as follows:

.. math::
   k_{eq} = \frac{ k_{cat1}k_{ip}k_{q} }{k_{cat2}k_{ia}} = \frac{ k_{cat1}k_{p}k_{iq} }{k_{cat2}k_{a}}

Second, the thermodynamic parameters in a network must not jointly imply that
it is possible to create or destroy energy simply by following a series of
reactions round in a loop. This implies that, at equilibrium, the net change in
Gibbs free energy due to the reactions in a loop should be exactly zero. In
mathematical notation:

.. math::
   \Sigma_{i\in loop}\Delta G_i = 0

Since there is a one-to-one relationship between :math:`k_eq` s and
:math:`DeltaG` s, this condition further constrains the feasible area of
thermodynamic parameter space for networks with loops.

How Maud ensures thermodynamic consistency
==========================================

In order to ensure that each enzyme's kinetic and thermodynamic parameters
agree, Maud ensures that one parameter per Haldane relationship is fixed based
on the values of the other parameters. For example, in the ordered unibi case
the :math:`k_{ip}` and :math:`k_{iq}` parameters are fixed as follows:

.. math::
   k_{ip} = \frac{k_{eq}k_{ia}k_{cat2}}{k_{q}k_{cat1}} \\
   k_{iq} = \frac{k_{eq}k_{cat2}k_{a}}{k_{cat1}k_{p}}

In order to avoid free energy loops, Maud generates :math:`k_{eq}` parameters
from :math:`\Delta G` parameters according to the following equation:

.. math::
   \mathbf{k_{eq}} = \exp(\frac{\Delta G}{-RT})

where R is the universal gas constant and T is the temperature in kelvin
(currently this is assumed to be 298). :math:`Delta G` parameters, in turn, are
generated as follows:

.. math::
   \Delta G = K\mathbf{b}


where :math:`\mathbf{b}` of auxiliary basis parameters whose length is the same
as the rank of the network's stoichiometric matrix and :math:`K =
Nullspace(Nullspace(S^{T})^{T})` is a matrix generated from the network's stoichiometric matrix
:math:`S` so as to ensure that :math:`\Delta G` sums to zero for loops.

In the case where the network has no loops, the width of `K` will be a diagonal
matrix and the basis parameters directly determine to the :math:`\Delta G` s. If
there are loops, there will be fewer basis parameters than :math:`\Delta G` s.

Information about marginal values of :math:`Delta G` - for example from *in
vitro* measurements - is represented directly as prior distributions on the
transformed :math:`Delta G` parameters. Since the transformation from basis
parameters to :math:`Delta G` s is linear and the posterior only needs to be
ascertained up to proportionality, there is no need for any adjustments to take
into account the effect of this transformation.
