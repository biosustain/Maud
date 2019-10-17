======================================
Thermodynamic parameters
======================================

Each enzyme in a reaction network model has a thermodynamic parameter called
:math:`k_{eq}`, which represents how much energy the corresponding chemical
reaction stores (:math:`k_{eq}` less than 1) or releases (:math:`k_{eq}`
greater than 1) from its environment. This document explains how the laws of
thermodynamics constrain these parameters and how Maud ensures that these
constraints are satisfied.

Thermodynamic constraints
=========================

Thermodynamic parameters are constrained in two main ways.

First, each thermodynamic parameter must agree with the kinetic parameters
describing its enzyme, according to the Haldane relationships governing the
enzyme's mechanism. For example, the Haldane relationships for an enzyme with
an ordered unibi mechanism are as follows:

.. math::
   k_{eq} = \frac{ k_{cat1}k_{ip}k_{q} }{k_{cat2}k_{ia}} = \frac{ k_{cat1}k_{p}k_{iq} }{k_{cat2}k_{a}}

Second, the thermodynamic parameters in a network must not jointly imply that
it is possible to create or destroy energy simply by following a series of
reactions round in a loop. In other words, that the product of all Keq
parameters in any loop must be exactly 1. This requirement can be represented
mathematically as the following Wegscheider condition:

.. math::
   N(S)^T\ln\mathbf{k_{eq}} = \mathbf{0}

In this equation :math:`N(S)` represents the right nullspace of the network's
stoichiometric matrix :math:`S`, i.e. the set of all flux vectors
:math:`\mathbf{v}` such that :math:`S\mathbf{v} = \mathbf{0}`.

How Maud ensures thermodynamic consistency
==========================================

In order to ensure that each enzyme's kinetic and thermodynamic parameters
agree, Maud ensures that one parameter per Haldane relationship is fixed based
on the values of the other parameters. For example, in the ordered unibi case
the :math:`k_{ip}` and :math:`k_{iq}` parameters are fixed as follows:

.. math::
   k_{ip} = \frac{k_{eq}k_{ia}k_{cat2}}{k_{q}k_{cat1}} \\
   k_{iq} = \frac{k_{eq}k_{cat2}k_{a}}{k_{cat1}k_{p}}

In order to avoid free energy loops, Maud sets :math:`k_{eq}` parameters
according to the following equation

.. math::
   \mathbf{k_{eq}} = \exp(L\mathbf{b})

where :math:`L` is a kernel of the set of possible solutions to the network's
Wegscheider equation, which is precomputed, and :math:`\mathbf{b}` is a vector
of random variables, whose length matches the width of :math:`L`. In the case
where the network has no loops, the width of `L` will be the same as the number
of reactions and the values of :math:`\mathbf{b}` will directly determine those
of :math:`\mathbf{k_{eq}}`. If there are loops, :math:`\mathbf{b}` will have
fewer elements than :math:`\mathbf{k_{eq}}`.


Jacobian adjustments
====================
It is important for Maud to be able to represent information about marginal
values of :math:`\mathbf{k_{eq}}` that may be available from *in vitro*
measurements. This is done by setting prior distributions on the transformed
:math:`\mathbf{k_{eq}}` values, and adjusting the distribution to take into
account the transformation.
