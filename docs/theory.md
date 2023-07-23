# Theory

This page explains the scientific and statistical theory behind Maud, and then
explains in broad terms how Maud implements that theory.

Look out for a paper setting this out in more detail in the near future!

(sec-the-big-picture)=
## The big picture

The overall aim of an analysis using Maud is to synthesise information from the
following sources:

* Quantitative measurements of a cell's metabolites, enzymes and fluxes.
* A mechanistic model that describes how each reaction in the metabolic network
  works in terms of some unknown parameters. Since Maud's model involves
  thermodynamics and kinetics, and because it sounds cool, we call it a
  "thermokinetic model".
* Information about the parameters in the model from sources other than the
  measurements, such as databases like {cite}`10.1093/nar/gkw952`.

These information sources are interesting on their own, but it is clear that
they are even better in combination. For example, there is not much gained by
accurately parameterising an enzyme's mechanism without knowing the values of
the parameters, and the more sources of parameter information the better!

By treating this problem as a case of Bayesian statistical inference, it is
possible to partition these information sources into separate sub-models:

* A measurement model describing what is learned from any set of measurements using
  conditional probability distributions with the form $p(y\mid\hat{y})$, where
  $y$ represents the observed and $\hat{y}$ the true values of the measured
  quantities.
* A generative model, encompassing our thermokinetic model, that describes how
  measurable quantities depend on parameters, i.e. a function
  $g:\Theta\rightarrow\hat{\mathfrak{y}}$ mapping any possible configuration of
  parameters $\theta$ to a set of measurable quantities $\hat{y}$.
* A prior model describing what is known about the unobserved parameters using
  an unconditional probability function with the form $p(\theta)$.

Given submodels in this form, and a set of measurements $y$ the posterior
distribution $p(\theta\mid y)$ synthesises the available information as we
want. Questions about the network can be formulated in terms of integrals over
the posterior distribution. For example, suppose that our thermokinetic model
says that the flux through a reaction is given by the function $flux:
\Theta\rightarrow\mathbb{R}$; our model's probability, given measurements $y$,
that the flux is greater than $x$ is then $\int_x^\inftyp(\theta\mid
y)d\theta$.

The rest of this page sets out how Maud fleshes out the three sub-models and
then briefly explains how it goes about calculating the required posterior
integrals.

(sec-thermokinetic-models-in-general)
## Thermokinetic models in general

A thermokinetic model describes the rates (usually referred to as 'fluxes') of
the chemical reactions in a metabolic network parametrically, in terms of the
concentrations of the network's metabolites, some kinetic parameters describing
the network's reactions, some thermodynamic parameters and some parameters
describing the network's boundary conditions. [](#eq-thermokinetic-general)
therefore captures the general form of a thermokinetic model.

$$
v = f_v(x, \theta)
$$(eq-thermokinetic-general)

In {eq}`eq-thermokinetic-general`, the term $x$ represents a vector of
metabolite concentrations, the term $\theta$ represents a parameter vector, and
$v$ is a vector of real numbers, one per reaction.

Thermokinetic models of cell metabolism are often contrasted with
constraint-based models. Constraint-based models do not describe metabolic
fluxes parametrically; instead they simply identify constraints that exclude
certain flux values.

Many thermokinetic models are possible, depending on how the modeller trades
off detail against other factors such as computational feasibility, model
identification and ease of curation.

One thing that almost all thermokinetic models have in common is the assumption
that the stoichiometric coefficient---that is, the number of molecules created
or destroyed---of every metabolite in every reaction is known and collected in
a stoichiometric matrix $S$ with a row for each metabolite and a column for
each reaction.

(sec-mauds-thermokinetic-model)=
## Maud's thermokinetic model

The thermokinetic model that we chose to implement in Maud decomposes into
factors as shown in [](#eq-flux-decomposition). The full model is presented
here for completeness and to illustrate the general nature of the equations
that we use, but the details can safely be skipped over.

$$
f_v(x, \theta) = Enzyme \cdot k_{cat} \cdot Reversibility \cdot Saturation \cdot Allostery \cdot Phosphorylation
$$(eq-flux-decomposition)

Each of the terms on the right hand side of [](#eq-flux-decomposition) is a
function of $x$ and $\theta$ representing a conceptually distinct aspect of how
an enzyme-catalysed reaction works.

The term $Enzyme$ is a vector of non-negative real numbers representing the
concentration of the enzyme catalysing each reaction.

The term $k_{cat}$ is a vector of non-negative real numbers representing the
amount of flux carried per unit of saturated enzyme.

The term $Reversibility$ is a vector of real numbers capturing the impact of
thermodynamics on the reaction's flux, as shown in
[](#eq-reversibility).
$$
\begin{aligned}
  Reversibility &= 1 - \exp(\frac{\Delta_{r}G' + RT \cdot S^T \ln(x)}{RT}) \\
  \Delta_{r}G' &= S^{T}\Delta_{f}G' + n F \psi
\end{aligned}
$$(eq-reversibility)

In [](#eq-reversibility), the term

* $T$ is the temperature in Kelvin (a number),
* $R$ is the gas constant (a number),
* $\Delta_rG'$ is a vector representing the Gibbs free energy change of each
reaction in standard conditions,
* $\Delta_fG'$ is a vector representing the standard condition Gibbs free
energy change of each metabolite's formation reaction, or in other words each
metabolite's 'formation energy'.
* $n$ is a vector representing the number of charges transported by each
reaction.
* $F$ is the Faraday constant (a number)
* $\psi$ is a vector representin each reaction's membrane potential (these
numbers only matter for reactions that transport non-zero charge)

Note that, for reactions with zero transported charge, the thermodynamic effect
on each reaction is derived from metabolite formation energies. This
formulation is helpful because, provided that all reactions' rates are
calculated from the same formation energies, they are guaranteed to be
thermodynamically consistent.

The term $n$ accounts for both the charge and the directionality. For
instance, a reaction that exports 2 protons to the extracellular space in the
forward direction would have -2 charge. If a negatively charged molecule like
acetate is exported in the forward direction, $n$ would be 1.

Notice how this way of modelling the effect of transported charge does not take
into account that the concentration gradient used by the transport is that of
the dissociated molecules. Thus, this expression is only correct for ions whose
concentration can be expressed in the model only in the charged form; e.g.,
protons, $K^+$, $Na^+$, $Cl^-$, etc.

The term $Saturation$ in equation [](#eq-flux-decomposition) is a vector of
non-negative real numbers representing, for each reaction, the fraction of
enzyme that is saturated, i.e. bound to one of the reaction's substrates. To
describe saturation we use [](#eq-saturation), which is taken from
{cite}`liebermeisterModularRateLaws2010`.

$$
\begin{aligned}
  Saturation_r &= a \cdot \text{free enzyme ratio} \\
             a &= \prod_{\text{s substrate}}\frac{x_s}{km_{rs}} \\
\text{free enzyme ratio} &= \begin{cases}
                  \prod_{\text{s sustrate}} (1 + \frac{x_s}{km_{rs}})^{S_sr}
                  + \sum_{\text{c inhibitor}}\frac{x_c}{ki_{rc}} & r\text{ irreversible} \\
                  -1
                  + \prod_{\text{s sustrate}} (1 + \frac{x_s}{km_{rs}})^{S_sr}
                  + \sum_{\text{c inhibitor}}\frac{x_c}{ki_{rc}}
                  + \prod_{\text{p product}} (1 + \frac{x_p}{km_{rp}})^{S_pr}  & r\text{ reversible}
                  \end{cases}
\end{aligned}
$$(eq-saturation)




The term $Allostery$ in [](#eq-flux-decomposition) is a vector of non-negative
numbers describing the effect of allosteric regulation on each reaction.
Allosteric regulation happens when binding to a certain molecule changes an
enzyme's shape in a way that changes its catalytic behaviour. We use
[](#eq-allostery), originally from {cite}`popovaGeneralizationModelMonod1975`,
to describe this phenomenon.

$$
\begin{aligned}
  Allostery_r &= \frac{1}{1 + tc_r \cdot (\text{free enzyme ratio}_r \cdot \frac{Qtense}{Qrelaxed})^{subunits}} \\
       Qtense &= 1 + \sum_{\text{i inhibitor}} \frac{x_i}{dc_{ri}} \\
     Qrelaxed &= 1 + \sum_{\text{a activator}} \frac{x_a}{dc_{ra}}
\end{aligned}
$$(eq-allostery)

Finally, the term $Phosphorylation$ in [](#eq-flux-decomposition) captures an
important effect whereby enzyme activity is altered due to a coupled process of
phosphorylation and dephosphorylation.

$$
\begin{aligned}
Phosphorylation_r &= (\frac{\beta}{\alpha + \beta})^{subunits} \\
          \alpha &= \sum_{\text{p phosphoylator}} kcat_{p} \cdot concp_p \\
           \beta &= \sum_{\text{p dephosphoylator}} kcat_{d} \cdot concd_d \\
\end{aligned}
$$(eq-phosphorylation)

(sec-steady-state-assumption)=
### Steady state assumption

As well as assuming that the fluxes in our metabolic network will behave as our
thermokinetic model expects, we also assume that the system was measured in a
steady state, so that every non-boundary metabolite's concentration was not
changing. This assumption is represented mathematically in [](#eq-steady).

$$
S\cdot f_v(x, \theta) = 0
$$(eq-steady)

The steady state assumption removes degrees of freedom equal to the rank of $S$
from our model, so that $v$ is constrained to lie in the right null space of S.
Usually, given $\theta$ it is possible can solve [](#eq-steady) to find a steady
state metabolite concentration $x_{steady}$.

(sec-summary-mauds-full-thermokinetic-model)=
### Summary: Maud's full thermokinetic model

Maud's full thermokinetic model starts with a parameter vector $\theta$. It then
calculates the steady state metabolite concentration vector $x_steady$ by
solving [](#eq-steady).

:::{note}

Note that in practice we do not solve [](#eq-steady) directly but instead use
ODE simulation - see section {ref}`sec-solving-the-steady-state-problem` for
details

:::

Metabolite concentrations are now available and can be compared with metabolite
measurements. To find flux values to compare with measurements we simply
calculate $v_{steady}=f_v(x_{steady}, \theta)$.

(sec-measurement-model)=
### Measurement model

For measurements of metabolite and enzyme concentrations we use independent
lognormal regression models:

$$
\begin{aligned}
  y_{x} &\sim LN(\ln(x_{steady}), \sigma_x) \\
  y_{enzyme} &\sim LN(\ln(Enzyme), \sigma_{enzyme})
\end{aligned}
$$(eq-conc-measurement-model-ln)

For measurements of steady state fluxes we use independent linear regression
models:

$$
\begin{equation}
  y_{v} \sim N(f_v(x_{steady}, \theta), \sigma_v)
\end{equation}
$$(eq-conc-measurement-model-n)

We ensure that flux measurements correspond to the flux modes of the measured
network, so that the same pathway is never measured twice. See chapter 10 of
{cite}`palssonSystemsBiologyConstraintbased2015`, entitled "the left null space" for a
discussion of this issue.

We assume that the measurement errors $\sigma_x$, $\sigma_{enzyme}$ and
$\sigma_v$ are known for each measurement.

These measurements $y_x$ and $y_{enzyme}$ are typically derived from
quantitative metabolomics and proteomics experiments. Our choice of measurement
model for these measurements is imperfect in at least these ways:

* The measurement error is not known, and is in fact very difficult to
  estimate.
* The measurements are not independent, as they are typically far more reliable
  as to differences---both between metabolites in the same experiment and
  between the same metabolite in different experiments---than they are as to
  absolute concentrations.

These problems also apply to flux measurements, but there is another issue:
flux meausrements are derived from k

(sec-prior-model)=
### Prior model

For kinetic parameters including $km$, $kcat$, $ki$ and parameters governing
regulation we use independent informative lognormal prior distributions based
on information gleaned from online databases, literature searches and
intuition.

For boundary conditions including unbalanced metabolite concentrations,
boundary fluxes and enzyme concentrations we use informative independent
lognormal or normal prior distributions depending on the natural sign
constraints of the variable (for example fluxes are not constrained to be
positive, so in this case we use independent normal prior distributions). We
sometimes use measurements to determine informative prior distributions for
boundary conditions.

For thermodynamic parameters---i.e. metabolite formation energies---we use an
informative multivariate normal distribution that is derived from equilibrium
constant measurements reported in the NIST TECRDB
{cite}`goldbergThermodynamicsEnzymecatalyzedReactions2004`. The distribution is
calculated using the component contribution method
{cite}`noorConsistentEstimationGibbs2013` as implemented in the software equilibrator
{cite}`beberEQuilibratorPlatformEstimation2021`. In future we would like to use our
own software to generate these priors.

(sec-implementation)=
## Implementation

In order to implement the statistical model described above, Maud performs
posterior sampling using adaptive Hamiltonian Monte Carlo as provided by
[Stan](https://mc-stan.org/).

To our knowledge this is the only viable approach. Realistic models have too
many parameters for many Bayesian computation approaches, including rejection
sampling, ABC, Metropolis-Hastings and Gibbs sampling, while a similar study
{cite}`st.johnBayesianInferenceMetabolic2018`, corroborated by our experience,
indicates that variational inference is unlikely to provide satisfactory
approximations in this case.

We believe that the non-linear, multi-parameter equations shown in
{ref}`sec-mauds-thermokinetic-model` create particular problems for MCMC
sampling because they induce complex correlations in the posterior
distribution. Consequently, although adaptive Hamiltonian Monte Carlo makes
Bayesian inference for realistic thermokinetic models possible, many leapfrog
steps are typically required per sample, as the sampler must traverse the
posterior distribution in small steps in order to avoid discretisation errors.

(sec-solving-the-steady-state-problem)=
### Solving the steady state problem

As mentioned in {ref}`sec-summary-mauds-full-thermokinetic-model`, Maud's
generative model calculates steady state metabolite concentrations given
parameter values by solving the steady state equation [](#eq-steady). Since
adaptive Hamiltonian Monte Carlo requires gradients of the posterior
distribution, it is also necessary to calculate sensitivities of the steady
state solution with respect to all parameters.

Our approach to this problem is to choose a starting concentration vector
$x_{0}$ and a simulation time $t$, then find $x_t$ using numerical ODE
integration. To verify whether $x_t$ is a steady state we then evaluate $S\cdot
f_v(x_t, \theta)$ and check if the result is sufficiently close to zero.

Stan provides an interface to the Sundials ODE solver CVODES, including
gradient calculations.

For the systems we have investigated, this method works better than solving the
steady state problem using an algebra solver.

(sec-analysis-of-results)=
### Analysis of results

After sampling is complete, Maud uses the Bayesian analysis library
[arviz](https://python.arviz.org/en/0.14.0/index.html) to transform the results
into an
[InferenceData](https://python.arviz.org/en/0.14.0/api/inference_data.html)
object and save it as a json file, along with a range of files containing debug
information. These files can be used to validate the computation and to draw
conclusions about the measured system.


```{bibliography}
```
