# 18. reverting from ajdoint to bdf solver

Date: 2021-11-31

## Status

In Review

## Context
To determine the steady state conditions of the model, we require an ODE solver.
We switched to the adjoint ODE solver for faster gradient evaluations,
however, recent use suggests that the solver is unstable and fails silently.

## Decision
Reverting to the ode_bdf_tol() solver implemented in Stan still solves
stiff ODE problems but slightly slower. This is overlooked as the current
adjoint solver is too unstable to use and fails on larger systems. The current
example files in the `tests/data/` folder work appropriately, however, larger
systems fail with the step size approaching 0.

An attempt was made to make the flux calculations in the ode more stable by
using built in functions and converting products to sums of logs etc.
This did not help the adjoint solver and hence this revertion was made.

Future versions of Stan can easily be tested by reverting to the previous
adjoint solver specification and the input files will still accept the tolerances
for the solver.

## Consequences
Slightly slower but more stable ODE evaluation.
By tuning the ODE to increase tolerances without causing failures the speed is
comparable.