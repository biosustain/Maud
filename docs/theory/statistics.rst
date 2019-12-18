==================
Statistical issues
==================

This document addresses statistical questions to do with Maud. Specifically, it
sets out what Maud does from a statistical point of view, what statistical
assumptions it makes and how these assumptions can be evaluated.

Bayesian inference for metabolic networks
=========================================

In general Bayesian inference offers two ways of expressing information
probabilistically: 'prior' probability distributions over possible values of
unknown parameters and measurement models or 'likelihood functions' that
noisily relate parameter values to measurements. Together these determine a
posterior distribution that synthesises the two sources of information and can
be used for statistical inference.

Maud's parameters are as follows:

- thermodynamic parameters (the formation energy for each compound in the
  network)
- kinetic parameters governing the behaviour of each reaction under fixed
  conditions
- concentrations of unbalanced metabolites (i.e. metabolites whose
  concentrations can be consumed or produced when the system is in a steady
  state)
- concentrations of enzymes

The prior distributions for formation energies come from an upstream analysis
involving in vitro measurements and inference about the composition of each
compound. Kinetic parameter priors are inferred from online databases using the
'respectful modelling' paradigm. Priors for concentrations of unbalanced
metabolites and enzyme concentrations are determined respectively by
metabolomics and proteomics measurements where these are available and by
reasonable physical limits otherwise.

The numbers that Maud treats as measurements are metabolomics data relating to
balanced metabolites and fluxomics data meausuring the steady-state rates of
reactions. Balanced metabolite concentration measurements are treated as draws
from a log-normal distribution with known standard deviation parameter and mean
parameter :math:`\log(true\ concentration)`. Flux measurements are similarly
treated as draws from a normal distribution with known standard deviation and
mean :math:`true flux`.

Overall Maud is fairly straightforward from a statistical point of view. The
main difficulties are (a) specifying how parameter values determine true steady
state concentrations and fluxes and (b) making the resulting problem
computationally tractable despite the need to solve a large system of
non-linear equations.

Statistical assumptions
=======================

It is somewhat unusual to use metabolomics data to specify prior distributions
directly rather than adding an extra measurement model. This choice was made
for the sake of simplicity.

Maud's prior distributions are meant to represent the range of plausible values
that its parameters can have. This assumption can be evaluated by using prior
predictive checks can be used to assess whether the joint prior distribution
implies any implausible measurements, or fails to accommodate any plausible
possibilities.

Maud's measurement model assumes that the accuracies of its measurements are
known in advance and follow the specified distributions. This assumption can be
evaluated by generating new simulated measurements from the posterior
predictive distribution and seeing whether they systematically differ from the
observed measurements.


Comparison with similar work
===========================

Maud's statistical approach is essentially the same as the one explained in
[1]. The main difference from that paper is that Maud attempts to represent the
kinetics of each reaction exactly using equations derived from consideration of
the elementary mechanism of the relevant enzymes. In contrast, st. John et al
use linear-logarithmic kinetic laws. The tradeoff is exactness on the one hand
against computational tractability on the other.

References
==========

[1] St. John, P., Strutz, J., Broadbelt, L. J., Tyo, K. E. J., & Bomble, Y. J. (2018). Bayesian inference of metabolic kinetics from genome-scale multiomics data. bioRxiv, http://dx.doi.org/10.1101/450163
