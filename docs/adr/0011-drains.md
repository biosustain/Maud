# 11. Drains as boundary condtitions

Date: 2020-10-19

## Status

Approved

## Context

Drains are an important aspect of any model, they are an essential boundary
condition for cases such as biomass drains. Curently rate laws are only specified using
the Generalised MWC format with the catalytic aspect using the modular rate law.

Introducing drains into Maud requires implementing this into the ODE system which
relies on enzymes as the defining feature, which are attributes of reactions in Maud.
To specify drains we create a new class independent of reactions, despite the fact that
they occur in the system of ODEs. This benefits post-processing as some techniques rely 
on varying enzyme concentrations, which we've established that drains do not have. 

Drains are considered as priors to ensure that prior predictive checks are informative
of the experimental condition.

## Decision

Drains will not be considered as a reaction class.


## Consequences

Processing data and setting fixed drains will be easy. But the current implementation
requires all drains to be specified.
