Installation
============

This document explains how to install Maud.

Install
=======

We recommend using a fresh virtual environment to install Maud. To make one and
then activate it, run the following commands:

.. code-block:: console

    sudo pip3.7 install virtualenv     # if virtualenv isn't installed already
    python3.7 -m virtualenv maud_venv  # choose any name you like!
    source maud_venv/bin/activate

To install Maud and its python dependencies to your new virtual environment,
run this command:

.. code-block:: console

    pip install -e git+ssh://git@github.com/biosustain/Maud.git#egg=maud

This step requires that you have ssh access to the Maud github repository
This is unavoidable while the repository remains private, but should be
achievable if you can see this readme.

Cmdstanpy depends on `cmdstan <https://github.com/stan-dev/cmdstan>`_, which
needs to be installed too. Fortunately, cmdstanpy comes with a command line
script that installs cmdstan, so this step is pretty simple:

.. code-block:: console

    install_cmdstan

