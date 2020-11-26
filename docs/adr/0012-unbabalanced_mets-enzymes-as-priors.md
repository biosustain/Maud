# 12. Unbalanced metabolites and enzyme concentrations as priors

Date: 2020-10-30

## Status

In Review

## Context

Unbalanced metabolites and enzyme concentrations are boundary conditions for our ODE
model. Experimental conditions are defined with respect to these values and drains, 
which are already defined as priors. Therefore, our prior knowledge about the
metabolic phenotype is defined as what is measured about the boundary conditions. This
decision aims to shift the measurements of the enzymes and unbalanced metabolites from
the likelihood evaluations to the prior information.

The benefit of treating priors in this way is that we define a prior on the phenotype
rather than all possible phenotypes. However, boundary conditions that are unmeasured
are still considered using weakly informative priors (read: within biologically relevant
boundaries).

## Decision

Unbalanced metabolites and enzyme concentrations can also be considered as prior distributions
rather than likelihood evaluations.


## Consequences

Prior predictive checks will be more informative. A seperate measurement model must be
constructed to identify the means and standard devations from the prior. Additionally,
large shifts from the prior distribution are an indication of measurement bias.
