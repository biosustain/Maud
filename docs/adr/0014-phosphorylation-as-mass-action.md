# 12. Implementing phosphorylation as linear mass action kinetics

Date: 2020-12-03

## Status

Draft

## Context

Phosphorylation is the process which deactivates metabolic enzymes.
The process is often conducted using kinases, which phosphorylate enzymes,
and phosphotases, which dephosphorylate the enzyme. This is a major metabolic
regulator and is essential for the accurate simulation of metabolic networks.


## Decision

To mininimise the number of parameters and reducing sampling time,
the linear rate law was selected (see [1] for a review). The linear mechanism 
is approximately correct when the following assumptions are satisfied:

* [Kinase] and [Phosphotase] << [Metabolic Enzyme],
* Rapid equilibrium binding between phosphorylation enzymes and metabolic enzyme,
* [Kinase] << dissociation constant for Kinase,
* [Phosphotase] << dissociation constant for Phosphotse,
* Competitive binding of phosphotases and kinases is negligable,
* The ATP/ADP ratio remains approximately constant,
* Phosphorylation and Dephosphorylation is an irreversible process.

Using these assumptions the steady state phosphorylated concentration is defined as:

`fraction_phosphorylated = alpha / (alpha + beta)`.

Where, alpha and beta correspond to the phosphorylation and dephosphorylation
rates, which are linear with respect to the kinase and phosphotase concentrations.

```
alpha = kcat * [Kinase]
beta = kcat * [Phosphatase]
```

The activity of the metabolic enzyme is proportional to the dephosphorylated
amount. To avoid situations where the kinase and phosphatase have the opposite
impact on the target enzymem, Maud will only refer to them as activating
and inhibiting enzymes.

## Consequences

Higher level protein responses are not currently supported:
* saturtion of target enzyme,
* competitive interactions between kinases and phosphatases
* allosteric regulation of kinases and phosphatases

## Literature

[1]	Salazar, C., & Höfer, T. (2006). Kinetic models of phosphorylation cycles: 
	A systematic approach using the rapid-equilibrium approximation for protein–protein 
	interactions. Biosystems, 83(2), 195–206. 
	https://doi.org/10.1016/j.biosystems.2005.05.015
