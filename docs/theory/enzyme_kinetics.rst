===============
Enzyme kinetics
===============

This document explains the assumptions about enzyme kinetics that Maud uses.

Modular rate law
================

This section outlines the assumptions made with the modular rate law, and includes a derivation of 2 substrate (A, B),
2 product (P, Q) random mechanism with competitive inhibitor I. It also highlights the general structure of the
modular rate law used in Maud. The modular rate law framework was taken from {Insert Modular rate law paper here},
and was adapted to suit our structure.

Assumptions
-----------

The assumptions used in the modular rate law are listed below:
    - The metabolite binding occurs in a random order,
    - Binding does not occur simultanesouly,
    - Substrates and products cannot bind at the same time,
    - Metabolite binding rates are much higher than the interconversion of substrate to product (Rapid Equilibrium),
    - Metabolite binding affinity is independent of order.


Example: 2 products and 2 substrate network
-------------------------------------------

For a random Bi-Bi network with the above assumptions, the rate will be the following:

.. math::
   v = E_t \frac{kcat_1 a' b' - kcat_2 p' q'}{(1 + a')(1 + b') + (1 + p')(1 + q') -1}

where for metabolite X the corresponding term is given by,

.. math::
   x' = \frac{X}{K_m^{x}}

A derivation is shown below. The rate is determined by the interconversion between substrate to product and using
elementary mass action kinetics is:

.. math::
   v = kcat_1 EAB - kcat_2 EPQ

Because of the rapid equilibrium assumption and dissociation constants are assumed 
for the Michaelis Menten constants EAB and EPQ can be written in terms of free
enzyme concentration and metabolite concentration. In this case:

.. math::
   EA = a' E_0
   EB = b' E_0
   EAB = a' EB, or,
   EAB = b' EA
   EAB = a' b' E_0

   EP = p' E_0
   EQ = q' E_0
   EPQ = p' EQ, or,
   EPQ = q' EP
   EPQ = p' q' E_0

After substituting the pseudo steady state enzyme concentration into the rate
equation it becomes:

.. math::
   v = E_0 (kcat_1 a' b' - kcat_2 p' q')

where the free enzyme ratio E_0 is defined by:

.. math::
   E_0 = \frac{1}{(1 + a')(1 + b') + (1 + p')(1 + q') -1}

Allostery
---------

Differing from the modular rate law defined in {Insert Modular rate law paper here},
allostery is considered using the generalised MWC form [see allostery link]. This 
requires the free enzyme amount - calculated above.