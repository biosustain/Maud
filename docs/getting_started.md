# Getting started

This page explains how to install Maud and use its command line interface.

## Installing Maud

Maud is compatible with Python versions 3.10 and above.

We recommend using a fresh virtual environment to install Maud. To make one and
then activate it, run the following commands:

```sh
python -m venv .venv --prompt=maud
source .venv/bin/activate
```

To install the latest version of Maud and its python dependencies to your new
virtual environment, run this command:


```sh
pip install maud-metabolic-models
```

Now Maud is installed!

### Installing cmdstan

Maud depends on the C++ application [cmdstan](
https://mc-stan.org/users/interfaces/cmdstan), which in turn requires a c++
toolchain.

To save users the hassle of installing one of those, we ship Maud with
pre-built binary files that make it just work on many computers. To see if your
computer is one of these, try running this command after installing Maud:

```sh
maud simulate
```
If some output appears, then you don't need to do anything!

If you see an error message, you will need to install cmdstan in order to use
Maud. See
[here](https://cmdstanpy.readthedocs.io/en/v1.1.0/installation.html#cmdstan-installation)
for a guide to how to do this using the Python package
[cmdstanpy](https://cmdstanpy.readthedocs.io/), which should be installed after
running `pip install maud-metabolic-models`.

## Using Maud

Maud is intended to be used by running commands starting with `maud` from the command line.

To see the available commands, type

```sh
maud --help
```
Maud's main commands take a path to an input folder as their only argument. See
the [inputting](#inputting) page for how to create input folders.

Maud provides a few pre-made inputs and the command `maud load-input <name>`
for loading them. You can see the available inputs by running `maud load-input
--help`. Try loading the `linear` input as with this command:


```sh
maud load-input linear
```

A typical workflow for using Maud, after creating an input folder, begins by
running the command `maud simulate path_to_my_input`. This will run
Maud's statistical model in simulation mode, generating a single snapshot of
the metabolic steady state corresponding to the initial parameter values. This
is useful for checking that the model and initial parameter values are sensible
before starting a sampling run.

Try doing this with the default input by running this command:

```sh
maud simulate linear
```

Maud should display some interesting information to your console. It should
also create a folder whose name starts with `maud_output_sim_linear`.
This contains all the information you might need for downstream analysis.

The next step is to run the command `maud sample` to trigger a sampling
run. You can do this as follows:

```sh
maud sample linear
```

When the sampler finishes, another folder should appear with a name beginning
`maud_output_linear`.

To make out of sample predictions, given a sampling run you can use the command
`maud predict` as follows, replacing `<your-timestamp>` with the
actual timestamp:

```sh
maud predict maud_output_linear<your-timestamp>
```

When this finishes the output folder should be augmented with prediction
information.
