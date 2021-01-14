# 13 SBML input

Date: 2021-14-01

# Context

Taking in an SBML model would make Maud a bit easier to pitch to people as it would remove the need to learn our custom syntax. I think it's possible to represent everything we need for a kinetic model: see [here](https://gist.github.com/teddygroves/287bd45c2cdc2bbb2defb4b47ff8fe53) for a quick example. There is still a little learning overhead as modifiers, drains and enzyme-catalysed reactions need to be appropriately annotated, and species names have to have the format `"{metabolite}_{compartment}"`. Also the user should be aware that Maud will ignore any attempts to specify rate equations. Nonetheless overall this seems to me like a better format for kinetic models.

However, SBML isn't a natural format for the quantitative parts of a Maud input, i.e. experiment results and priors. There is a fairly recent extension to SBML called [distributions](http://co.mbine.org/specifications/sbml.level-3.version-1.distrib.version-1.release-1) that is meant to handle probability distributions, but I don't think it's a good fit for our purposes as the specification is quite restrictive (for example you can't define a lognormal distribution using quantiles and multivariate distributions aren't allowed) and the focus seems to be on creating a setup for prior sampling rather than simply recording measurements and priors.

I think the best solution is to allow SBML specification of kinetic models, but to keep these in different files from the quantitative information. We have been moving in this general direction for a while - this was the long term reason for specifying the kinetic model, priors and experiments separately - and now seems like a good time to cut the cord.

In more detail, this would mean changing how `maud sample` works so that it requires separate kinetic model and quantitative information arguments, e.g

```
maud sample kinetic_model.sbml quantitative_information.toml
```

or perhaps

```
maud sample kinetic_model.sbml --experiments experiments.toml --priors priors.json
```

# Decision
I think I prefer the first approach as in practice experiments and priors are quite coupled so it makes sense to define them in the same file.

## Consequences
The main work involved in implementing this proposal is to update the `cli` and `sampling`modules in order to accommodate multiple files and also to update the `io` module so it can parse sbml format kinetic models. I think these should both be fairly straightforward.
