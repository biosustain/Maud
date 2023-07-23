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

Maud is an application that fits Bayesian statistical models of
metabolic networks using `Python <https://www.python.org/>`_ and `Stan
<https://mc-stan.org>`_.

Maud aims to take into account allosteric effects, ensure that the laws of
thermodynamics are obeyed and to synthesise information from both steady state
experiments and the existing literature.

Installation
============
First create a fresh Python virtual environment and then activate it:

.. code-block:: console

    python -m venv .venv --prompt=maud
    source .venv/bin/activate

To install Maud and its python dependencies to your new virtual environment, run
this command:

.. code-block:: console

    pip install maud-metabolic-models

Cmdstanpy depends on `cmdstan <https://github.com/stan-dev/cmdstan>`_, 
which in turn requires a c++ toolchain. On some computers you will have to 
install these in order to use Maud. You will hit an error at the next step if 
this applies to your computer. Luckily cmdstanpy comes with commands that
can do the necessary installing for you. On windows the c++ toolchain can be installed with 
the following powershell commands:


Usage
=====
Maud is used from the command line. To see all the available commands try running 

.. code-block:: console

    maud --help

Copyright
=========

* Copyright (c) 2023, Novo Nordisk Foundation Center for Biosustainability, Technical University of Denmark.
* Free software distributed under the `GNU General Public License 3.0 <https://www.gnu.org/licenses/>`_
