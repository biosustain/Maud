# 17. conducting in and out of sample simulations of the posterior distribution

Date: 2021-01-22

## Status

Draft

## Context
Currently, simulations are limited to trained data, with no functionality 
to include out-of-sample simulations using the posterior parameters. The goal
is to take posterior samples and simulate both the trained data and additional experiments.
This will allow validation of the posterior samples, but is more flexible than
definig the validation data at the time of sampling. To accomplish this we will
introduce a new stan program "simulate_ensemble.stan", which, will take posterior
samples and priors on boundary conditions as input. We can then sample from the
boundary conditions on "unbalanced_metabolites", "drains", and "enzyme_concentrations" 
for the out-of-sample dataset.

## Decision
To include a seperate function that is directed to a folder containing draws
from the same kinetic model, which, may or may not contain the same experimental conditions.
We will then parse this to a stan model performing simulations in the generated_quantities
block, conducted using the "fixed_params=True" argument for the cmdstanpy.CmdStanModel.sample
object. Redefining the variables in the generated_quantities block will make processing
data much easier. The experiments will need to be split according to if they were sampled
or not; because, this will determined if posterior draws or prior samples are used for
simulations.