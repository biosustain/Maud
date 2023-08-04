Welcome to Maud's documentation!
================================

Maud is a tool for fitting Bayesian statistical models of metabolic networks.

Maud's distinguishing features include:

* Scientifically accurate representation of phenomena including enzyme kinetics,
  allosteric regulation, competitive inhibition, phosphorylation, knockouts and
  transported charges.
* Guaranteed consistency with thermodynamic and steady state constraints.
* A statistical model allowing inference consistent with both information from
  experiments and pre-experimental background information.
* Prediction of steady state concentrations and fluxes given out-of-sample
  boundary conditions.

More practically speaking, Maud is a Python application that you can use from
the command line. Maud provides the command `maud sample`, which loads data
from an input folder containing [`toml`](https://toml.io/en/) files
representing a metabolic model, experimental measurements and background
information, validates this input and passes it to a statistical model
specified in using the probabilistic programming language [Stan](
https://mc-stan.org/). When the statistical computation is finished, Maud
converts the results into a convenient format and saves them in an output
folder.

To find out how to install and use Maud, check out the
[](#getting_started) page.

For a guide to creating inputs for Maud, see the [](#inputting) page.

For what to do when things inevitably go wrong, read the
[](#troubleshooting) page.

For detailed discussion of the scientific assumptions implicit in Maud's
statisical model, there is the [](#theory) page.

For reference information about Maud's command line interface, see the
[](#cli) page.

To find out how to contribute to Maud, read the [](#contributing)
page.

```{toctree}
:maxdepth: 0
:numbered:

getting_started.md
inputting.md
troubleshooting.md
theory.md
contributing.md
cli.md
```
