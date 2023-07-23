====================
Maud's Documentation
====================

This note explains how to build Maud's documentation locally and how to
contribute to it.

Building the docs locally
=========================

In order to build the documentation locally you will need to have `sphinx
<http://www.sphinx-doc.org/en/master/>`_ installed in your python
environment. If you don't, you should be able to install it with the following
command line incantation:

.. code:: bash

    pip install -U Sphinx

Sphinx turns the source files in the `docs
<https://github.com/biosustain/Maud/tree/master/docs>`_ directory into nicely
formatted and structured output files in the directory :code:`docs/_build`. You
can check all the available options by running :code:`make` from the docs
directory.

In order to create a set of html pages, run

.. code:: bash

    make html

This will create an html version of the documentation, which you can read by
opening the file :code:`docs/_build/html/index.html` in your browser.

Alternatively, to create a pdf version of the documentation, run

.. code:: bash

    make latexpdf

The output should appear at :code:`docs/_build/latex/maud.pdf`.


Contributing to the documentation
=================================

To add new documentation, either edit an existing :code:`.rst` file in the
docs directory, or add a new one and edit the file `docs/index.rst
<https://github.com/biosustain/Maud/blob/master/docs/index.rst>`_ so that the
new file is included in the :code:`toctree` list.

To check that the new contribution looks like you expect, you can build the
docs locally using the commands :code:`make html` or :code:`make latexpdf` as
described above.
