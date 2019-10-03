=============================
Maud
=============================

.. image:: https://img.shields.io/pypi/l/maud.svg
   :target: https://www.apache.org/licenses/LICENSE-2.0
   :alt: Apache Software License Version 2.0

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/ambv/black
   :alt: Black


Maud is an application that fits Bayesian statistical models of metabolic
networks. It aims to take into account allosteric effects, ensure that the laws
of thermodynamics are obeyed and to synthesise information from both steady
state experiments and the existing literature.

Install
=======
First create a fresh Python 3.7 virtual environment and then activate it:

.. code-block:: console
    sudo pip3.7 install virtualenv     # if virtualenv isn't installed already
    python3.7 -m virtualenv maud_venv  # choose any name you like!
    source maud_venv/bin/activate

Run this command to install Maud and its python dependencies:

.. code-block:: console
    pip install -e git+ssh://git@github.com/biosustain/Maud.git#egg=maud

This step requires that you have ssh access to the Maud github repository
This is unavoidable while the repository remains private, but should be
achievable if you can see this readme.

Cmdstanpy depends on `cmdstan <https://github.com/stan-dev/cmdstan>`_, which needs to be installed too. Fortunately,
cmdstanpy comes with a command line script that installs cmdstan, so this step
is pretty simple:

.. code-block:: console
    install_cmdstan


Usage
=====
To run the simple linear model, use the following command:

.. code-block:: console
    maud sample

This will use the data file at `data/in/linar.toml` to create a Stan program at
`maud/stan_code/autogen/inference_model_linear.stan`, compile it into a
C++ Stan model, draw samples from the resulting posterior and store them in
files starting with `data/out/model_output_linear`.

The `sample` command can be configured in a few ways - to check out all the
options try running

.. code-block:: console
    maud sample --help


Copyright
=========

* Copyright (c) 2019, Novo Nordisk Foundation Center for Biosustainability, Technical University of Denmark.
* Free software distributed under the `GNU General Public License 3.0 <https://www.gnu.org/licenses/>`_
