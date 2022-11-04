.. Maud documentation master file, created by
   sphinx-quickstart on Fri Oct  4 15:14:47 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Maud's documentation!
================================

Maud is a tool for fitting Bayesian statistical models of metabolic networks.

Maud's distinguishing features include:

- Scientifically accurate representation of phenomena including enzyme kinetics,
  allosteric regulation, competitive inhibition, phosphorylation, knockouts and
  transported charges.
- Guaranteed consistency with thermodynamic and steady state constraints.
- A statistical model allowing inference consistent with both information from
  experiments and pre-experimental background information.
- Prediction of steady state concentrations and fluxes given unseen boundary
  conditions.

More practically speaking, Maud is a command line application that uses Stan to
specify and fit a statistical model and Python to interface between Stan and
humans.

.. toctree::
   :maxdepth: 1
   :caption: How to use Maud:

   usage/installation
   usage/inputting
   usage/post_installation_usage
   usage/contributing

.. toctree::
   :maxdepth: 1
   :caption: Theoretical background:

   theory/enzyme_kinetics
   theory/kinetic_model
   theory/statistics
   theory/thermodynamics
   theory/drains
   theory/transported_charges
   theory/papers

..
   .. toctree::
      :maxdepth: 1
      :caption: Case studies:

      case_studies/ecoli


.. toctree::
   :maxdepth: 1
   :caption: How Maud works:

   how_maud_works/data_model
   how_maud_works/computation



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
