# 17. conducting out of sampling 

Date: 2021-08-31

## Status

partially implemented

## Context
After sampling from the posterior distribution users may want the option
of validating their results against experimental data or to predict 
cellular behaviour using a trained kinetic model.

## Decision
Posterior predictive samples will be implemented using the posterior
draws for kinetic and thermodynamic parameters. The boundary conditions
from the predicted experiments will be sampled from their marginal distributions.
The number of samples will be limited to the number of samples from the posterior
distribution.

In order to define what is part of the training set and prediction set a new
file will be introduced called experiment_metadata. A toml file where the training
prediction split will be defined under the headers ["training"], and ["prediction"]
respectively.

A new stan file will be made where there is no model block (minimal if required).
In the generated_quantities (gqs), we will iterate through the draws from a previously
generated set of csvs. the drains, conc_enzymes, and conc_unbalanced_metabolites 
will be sampled from their priors using the functions normal_rng() and lognormal_rng()
for their respective generative distributions.

log_probabilities should be calculated if measurements are also included for the
prediction experiments.

## Consequences
The inclusion of a new input file would mean previous inputs would be broken.
However, these fixes are relatively small as it's only a list.

On the other hand you can now perform validation and prediction - which
expand the Maud toolset dramatically.