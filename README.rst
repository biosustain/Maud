====
Maud
====

.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0
   :alt: GNU General Public License 3.0

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/ambv/black
   :alt: Black

.. image:: https://img.shields.io/badge/Contributor%20Covenant-v1.4%20adopted-ff69b4.svg
   :target: https://www.contributor-covenant.org/
   :alt: Contributor Covenant Version 1.4

.. image:: https://readthedocs.org/projects/maud-metabolic-models/badge/?version=latest
   :target: https://maud-metabolic-models.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

Maud is a work-in-progress application that fits Bayesian statistical models of
metabolic networks using `Python <https://www.python.org/>`_ and `Stan
<https://mc-stan.org>`_.

Maud aims to take into account allosteric effects, ensure that the laws of
thermodynamics are obeyed and to synthesise information from both steady state
experiments and the existing literature.

Install
=======
First create a fresh Python 3.7 virtual environment and then activate it:

.. code-block:: console

    sudo pip3.7 install virtualenv     # if virtualenv isn't installed already
    python3.7 -m virtualenv maud_venv  # choose any name you like!
    source maud_venv/bin/activate

To install Maud and its python dependencies to your new virtual environment, run
this command:

.. code-block:: console

    pip install maud-metabolic-models

Cmdstanpy depends on `cmdstan <https://github.com/stan-dev/cmdstan>`_, 
which in turn requires a c++ toolchain. Fortunately, cmdstanpy comes with
commands that can install these for you. On windows the necessary dependencies 
can be installed with the following powershell commands:

.. code-block:: console

    python -m cmdstanpy.install_cxx_toolchain
    python -m cmdstanpy.install_cmdstan --compiler

On macos or Linux, you must install the c++ requirements manually 
(see [here](https://cmdstanpy.readthedocs.io/en/v0.9.67/installation.html#install-cmdstan)) for instuctions. 
Cmdstan can then be installed using this shell command:

.. code-block:: console

    install_cmdstan

Usage
=====
To run the simple linear model, use the following command:

.. code-block:: console

    maud sample

This will compile the Stan program at `src/maud/inference_model.stan
<https://github.com/biosustain/Maud/blob/master/src/maud/inference_model.stan>`_, 
then run the resulting binary file using the data at `tests/data/linar.toml
<https://github.com/biosustain/Maud/blob/master/tests/data/linear.toml>`_, storing
the results in csv files starting with
:code:`model_output_linear`.

The `sample` command can be configured in a few ways - to check out all the
options try running

.. code-block:: console

    maud sample --help


Copyright
=========

* Copyright (c) 2019, Novo Nordisk Foundation Center for Biosustainability, Technical University of Denmark.
* Free software distributed under the `GNU General Public License 3.0 <https://www.gnu.org/licenses/>`_
