# 18. reverting from ajdoint to bdf solver

Date: 2022-03-29

## Status

In Progress

## Context

Rational engineering of strains is a difficult problem due to the non-linearity
of metabolism. Metabolic Control Analysis (MCA) assists scientists by providing
them with flux control coefficient. flux control coefficients state how changing
a parameter will influence the flux for any given reference state.


## Decision
Calculate the flux control coefficients and metabolite control coefficients for
every condition that is sampled.


## Consequences

Good Consequences:
Users will have access to MCA for every simulation that they run.


Bad Consequences:
Larger file sizes that scales by:
(n_experiments*n_fluxes*n_enzymes)+(n_experiments*n_metabolites*n_fluxes).
This also requires Stan 2.29+ as the implementation relies on complext numbers.
