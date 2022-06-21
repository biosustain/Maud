============
Installation
============

This page explains how to install Maud.

We recommend using a fresh virtual environment to install Maud. To make one and
then activate it, run the following commands:

.. code-block:: console

    python -m venv maud_venv  # choose any name you like!
    source maud_venv/bin/activate

To install the latest Maud and its python dependencies to your new virtual
environment from the latest master branch, run this command:

.. code-block:: console

    pip install https://github.com/biosustain/Maud/archive/master.zip

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

