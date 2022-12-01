============
Installation
============

This page explains how to install Maud.

Maud is compatible with Python versions 3.9 and above

We recommend using a fresh virtual environment to install Maud. To make one and
then activate it, run the following commands:

.. code-block:: console

    python -m venv .venv --prompt=maud
    source .venv/bin/activate

To install the latest Maud and its python dependencies to your new virtual
environment, run this command:

.. code-block:: console

    pip install maud-metabolic-models

Cmdstanpy depends on `cmdstan <https://github.com/stan-dev/cmdstan>`, which
in turn requires a c++ toolchain. Fortunately, cmdstanpy comes with commands 
that can install these for you. On windows the necessary dependencies can be 
installed with the following powershell commands:

python -m cmdstanpy.install_cxx_toolchain
python -m cmdstanpy.install_cmdstan --compiler

On macos or Linux, you must install the c++ requirements manually 
(see [here](https://cmdstanpy.readthedocs.io/en/v0.9.67/installation.html#install-cmdstan))
for instuctions. Cmdstan can then be installed using this shell command:

.. code-block:: console

    install_cmdstan

