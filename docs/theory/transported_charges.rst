==============================
Transport of charged molecules
==============================

This document explains the assumptions about the transport of charged molecules
in Maud.

Transport reactions
===================

The driving force of membrane transport is affected by the membrane potential
in the form:

.. math::
	\Delta G_{transport} = \Delta G_r^0 + RT S log(conc) + n F \psi

where :math:`\Delta G_r^0` is the standard free reaction energy from formation
energies; :math:`RT` is the gas constant times the temperature; :math:`S` is
the stoichiometric matrix; :math:`F` is the Faraday constant; :math:`\psi` is
the membrane potential; and :math:`n` is the number of charges being exchanged.
Notice how if there is no charge transport, :math:`n` is 0 and the driving
force expression becomes that of a normal driving force.

:math:`n` accounts for both the charge and the directionality. For instance, a
reaction that exports 2 protons to the extracellular space in the forward
direction would have -2 charge. If a negatively charged molecule like acetate
is exported in the forward direction, :math:`n` would be 1.

Notice how this does not take into account that the concentration gradient used
by the transport is that of the dissociated molecules. Thus, as of now, this
expression is only correct for ions whose concentration can be/is expressed in
the model only in the charged form; e.g., protons, :math:`K^+`, :math:`Na^+`,
:math:`Cl^-`, etc.

Implementation
==============

A prior of name :code:`psi` is used in the above equation to account for the membrane
potential. This prior is usually negative, around -0.95V and is tied to a
particular experiment. The value for the charge is a structural parameter of 
the reaction and the rest (stoichiometry, concentrations) work as for any other
reaction.

The number of transported charges :math:`n` can be explicited in the kinetic model
TOML file for each reaction as :code:`transported_charge` (0 by default). For
instance:

.. code:: toml

   [[reaction]]
   id = "ATPSm"
   name = "ATPase c0"
   mechanism = "reversible_modular_rate_law"
   water_stoichiometry = 1.0
   transported_charge = 3.3
   stoichiometry = { adp_c = -1.0, h_e = -3.3, h_c = 3.3, pi_c = -1.0, atp_c = 1.0 }
