====================
Contributing to Maud
====================

This document explains how you can contribute to Maud.


Reporting bugs and making feature requests
==========================================

If you find a bug in Maud or would like to request a new feature, consider
adding to the `issues list <https://github.com/biosustain/Maud/issues>`_. Bugs
should get the red :code:`bug` label and features the turquoise
:code:`enhancement` label.


Contributing code changes
=========================

If you want to add a new feature or fix a bug yourself, the first step is to
clone the code. If you have ssh access to the `Maud github repository
<https://github.com/biosustain/Maud>`_, you can clone it by running

.. code:: bash

          git clone git@github.com:biosustain/Maud.git

Next, check out a new branch

.. code:: bash

          git checkout -b descriptive_branch_name

When you are happy with your changes, commit them to your new branch and then
push it to github

.. code:: bash

          git commit -m "Short description of the changes I've made"
          git push

Finally, use the online github interface to make a new pull request from your
branch, complete with a longer description.

At least one approving review is required before a pull request can be
merged. You can make your reviewers' lives a lot easier by adding as much
detailed commentary as possible, referring to existing issues where appropriate
and keeping your pull requests small and logically distinct.


Running tests locally
=====================

You will need to have the following python packages installed in order to run
Maud's tests locally

* `tox <https://tox.readthedocs.io/en/latest/>`_
* `black <https://github.com/psf/black>`_
* `flake8 <http://flake8.pycqa.org/en/latest/>`_
* `isort <https://github.com/timothycrosley/isort>`_

Follow the links for installation instructions - they are mostly

.. code:: bash

          pip install <package_name>

To run the tests locally, run the command

.. code:: bash

          tox

To painlessly ensure that your code complies with black and isort's guidelines,
you can run this reformatting command:

.. code:: bash

          make qa

Releasing versions
==================

Maud uses `versioneer<https://github.com/warner/python-versioneer>`_ to handle
releases.

To release a new version of Maud, first make a pull request updating the
changelog with the new version number and explaining what has changed since the
previous version.

Once the changes are merged into the :code:`origin/master` branch, add a tag
with the new version number to your local :code:`master` branch:

.. code:: bash

          git tag 0.2.1

Now push the new tag to github:
 
.. code:: bash

          git push --tags
