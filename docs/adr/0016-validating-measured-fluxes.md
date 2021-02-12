# 16. checking that appropriate fluxes are measured

Date: 2021-01-22

## Status

In Review

## Context
Defining measurements for independent fluxes isn't always clear.
This can occur when you measure more fluxes than there are degrees
of freedom in a network.

An example would be this simplified network:

	A -> B -> C

where reaction 1 and reaction 2 are dependent, implying that
no additional information is achieved by including both.

Another issue is knowing when you do not have enough fluxes
measured, resulting in an underdetermined system. Due to the
Bayesian implementation of Maud, these systems are still theoretically
resolvable. However, supplementing as much information
as possible will likely be beneficial.

## Decision
Identifying underdetermined systems is acomplished by first calculating
the null space of the matrix. This gives the number of degrees of freedom
of the system as well. Then we calculate the reduced row echelon form of
the transpose of the null space. The resulting matrix represents the
independent flux pathways through the network as rows. If you take the
measured subset of reactions and there is a row containing no non-zero
entries then the system is not fully described using the current measurements.

Determining if the system is overspecified is achieved by comparing the
number of measurements to the degrees of freedom. If the number of measurements
is larger than the degrees of freedom then the system is overdetermined.

It is possible to both have an underdetermined system which is overspecified
by having multiple measurements on dependent paths. It is also possible to
recieve the warning that the system is overspecified by independent measurements.
For instance, a linear pathway where the influx and efflux are both measured.
This is still valid as they are independent measurements.

## Consequences
There are no direct consequences to the user, this should not effect Maud sampling in any
way. The output of this feature may warn the user about possible biases present
in the experiment definition. It is however completely up to the user to determine
if these warnings apply to their notwork or not.