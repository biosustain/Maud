===============
Enzyme kinetics
===============

This document explains the assumptions about enzyme kinetics that Maud uses.

Modular rate law
================

This section outlines the assumptions made with the modular rate law, and includes a derivation of 2 substrate (A, B),
2 product (P, Q) random mechanism with competitive inhibitor I. It also highlights the general structure of the
modular rate law used in Maud. The modular rate law framework was taken from [1],
and was adapted to suit our structure.

Assumptions
-----------

The assumptions used in the modular rate law are listed below:
    - the metabolite binding occurs in a random order,
    - binding does not occur simultanesouly,
    - substrates and products cannot bind at the same time,
    - metabolite binding rates are much higher than the interconversion of substrate to product (rapid equilibrium assumption),
    - metabolite binding affinity is independent of order.


Example: 2 products and 2 substrate network
-------------------------------------------

.. figure:: random-bibi.png

    A random mechanism with 2 products and 2 substrates with a slow conversion step. All of the reactant
    binding/release steps are in rapid equilibrium.

For a random Bi-Bi network with the above assumptions, the rate will be the following:

.. math::
   v = E_t \frac{kcat_1 a' b' - kcat_2 p' q'}{(1 + a')(1 + b') + (1 + p')(1 + q') -1}

where for metabolite X the corresponding term is given by,

.. math::
   x' &= \frac{X}{K_m^{x}} \\
   K_m^{x} &= \frac{[X] \bullet [E_{i-1}]}{[E_i]}

where :math:`E_{i-1}` is the enzyme state not bound to metabolite X, and,
:math:`E_i` is the enzyme state bound to metabolite X

A derivation is shown below. The rate is determined by the interconversion between substrate to product and using
elementary mass action kinetics is:

.. math::
   v = kcat_1 EAB - kcat_2 EPQ

Because of the rapid equilibrium assumption the dissociation constants are assumed 
for the Michaelis Menten constants. EAB and EPQ can be written in terms of free
enzyme concentration and total enzyme concentration and metabolite concentrations. In this case:

.. math::
   EA &= a' E_0  \\
   EB &= b' E_0  \\
   EAB &= a' EB = b' EA = a' b' E_0 \\
   EP &= p' E_0  \\
   EQ &= q' E_0  \\
   EPQ &= p' EQ = q' EP = p' q' E_0 \\
   E_0 &= E_t - \sum_{i} E_i \\
   E_0 &= E_t - E_0 (a' + b' + a' b' + p' + q' + p' q') \\
   E_0 / E_t &= \Theta \\
   \Theta &= \frac{1}{1 + a' + b' + a' b' + p' + q' + p' q'}
   

After substituting the enzyme concentrations into the rate
equation it becomes:

.. math::
   v = E_0 (kcat_1 a' b' - kcat_2 p' q')

where the free enzyme amount :math:`E_0` is defined by:

.. math::
   E_0 = \frac{E_t}{(1 + a')(1 + b') + (1 + p')(1 + q') -1}

Competitive inhibition
----------------------
In the following case we will consider competitive inhibition where an inhibitor
selectively binds to the free enzyme, preventing binding from either substrate or
product.

.. figure:: random-bibi-competitive.png

    A random mechanism with 2 products and 2 substrates with a slow conversion step.
    all metabolites including the inhibitor are in rapid equilibrium with the enzyme
    states.

As described in [1], competitive inhibition is accounted for in the denominator
term of the rate equation. It's easy to see how this occurs when you look at the free
enzyme concentration:

.. math::
   EI = i' E_0

therefore,

.. math::
    \Theta = \frac{1}{1 + a' + b' + a' b' + p' + q' + p' q' + i'}

which can then be substituted into the original rate equation with the form:

.. math::
   v = E_t \frac{kcat_1 a' b' - kcat_2 p' q'}{(1 + a')(1 + b') + (1 + p')(1 + q') + i' -1}

Allostery
---------

Differing from the modular rate law defined in [1],
allostery is considered using the generalised MWC form [see allostery link]. This 
requires the free enzyme amount - calculated above.

References
==========
[1] Liebermeister, W., Uhlendorf, J. & Klipp, E. Modular rate laws for enzymatic reactions: 
thermodynamics, elasticities and implementation. Bioinformatics 26, 1528–1534 (2010).