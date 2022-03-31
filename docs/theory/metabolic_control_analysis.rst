==========================
Metabolic Control Analysis
==========================

Introduction to Metabolic Control Analysis
==========================================

A large portion of kinetic model development has the goal of improving
metabolism. Unfortunately, users have to deal with the non-linearities that
these models present. A consequence of this is that logical behaviour in the
local sense when looking at a single enzyme, will not translate the the global
response of the system. Metabolic control analysis is a tool that gives a greater
understanding of metabolic pathways as a whole and the response of the system
to changes in parameters. The implementation in Maud will only focus on
flux control coefficients with chagnes in enzymes, and metabolite control
coefficients.

Definition
==========

To understand the impact on the network we can take the derivative at the
steady state defined by :math:`Nv=0`. Resulting in:

.. math::
    S\frac{\partial v}{\partial X} \frac{dX}{dp} + S\frac{\partial v}{\partial p} = 0

This equation serves as the basis for calculating the Concentration control coefficient.
By rearranging this we can isolate :math:`\frac{dX}{dp}:

.. math::
    \frac{dX}{dp} = -[S \frac{\partial v}{\partial X}]^{-1} S\frac{\partial v}{\partial p}

where the concentration control coefficient is defined as:

.. math::
    C^X = -[S \frac{\partial v}{\partial X}]^{-1} S

and the flux control coefficient is:

.. math::
    C^v = 1 + \frac{\partial v}{\partial X}C^X


Limitations
===========
We do not deal with conservation relationships currently, which can be described
as a group of metabolites summing to a particular value. If we were to implement
this feature, the MCA should be altered to include the link matrix that is
omitted as it is the identity matrix.
