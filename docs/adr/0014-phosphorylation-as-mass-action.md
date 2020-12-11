# 12. Implementing phosphorylation as linear mass action kinetics

Date: 2020-12-03

## Status

Draft

## Context

Phosphorylation is the process which deactivates metabolic enzymes,
the process is often conducted using kinases, which phosphorylate enzymes,
and phosphotases, which dephosphorylate the enzyme. This is a major metabolic
regulator and is essential for the accurate simulation of metabolic networks.


## Decision

There are many potential regulatory mechanisms to implement, as discussed in
https://doi.org/10.1016/j.biosystems.2005.05.015. with the following assumptions
the selected rate law is the linear rate law:

* [Kinase] and [Phosphotase] << [Metabolic Enzyme],
* Rapid equilibrium binding between phosphorylation enzymes and metabolic enzyme,
* [Kinase] << dissociation constant for Kinase,
* [Phosphotase] << dissociation constant for Phosphotse,
* Competitive binding of phosphotases and kinases is negligable.
* The ATP/ADP ratio
* Phosphorylation and Dephosphorylation is an irreversible process.

Using these assumptions the steady state phosphorylated concentration is defined as:

`fraction_phosphorylated = alpha / (alpha + beta)`.

Where, alpha and beta correspond to the phosphorylation and dephosphorylation
rates, which are linear with respect to the kinase and phosphotase concentrations.

```
alpha = kcat * [Kinase]
beta = kcat * [Phosphotase]
```

The activity of the metabolic enzyme is proportional to the dephosphorylated
amount. In situations where the phosphorylated enzyme is active, for the sake
of simplicity the kinases and phosphotases should be swapped in Maud. 
Any measurements of phosphorylated fraction should also be inverted using the
calculation `phosphorylated fraction = 1 - dephosphorylated fraction`.

## Consequences
