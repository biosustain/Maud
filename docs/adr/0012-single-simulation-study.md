# 12. Script implementing a single simulation study

Date: 2020-11-13

## Status

Draft

## Context

We need a way to run simulation studies in order to verify that our models work
and to answer questions about how informative particular measurements and
priors are.


## Decision

We will initially make a python script that implements a single simulation
study, targeting the model at `ecoli_small_experiments.toml`. Based on how
that goes, we can then think about how to make a more general system for
running arbitrary simulation studies.

It's better to creep towards a simulation study system like this rather than
jump straight to the general tool because it's quite a big job, and we don't
really know exactly what the general tool will look like yet.

## Consequences

We'll probably need to delete the script and repeat some work when making the
more general simulation study system. 
