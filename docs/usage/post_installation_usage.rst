==========================
Generating posterior draws
==========================

Maud's main purpose is to generate an ensemble of samples that collectively
represent the information contained in some input data. 

To generate samples for input data at :code:`path/to/myinput`, run the
following command:


.. code:: console

          maud sample path/to/myinput


This will trigger maud to generate an output folder in the directory where the
command is run. To save the output folder somewhere else, use the option
:code:`output_dir` as follows:

.. code:: console

    maud sample path/to/myinput --output_dir path/to/desired/output/dir

Here is the full documentation for the :code:`sample` command, listing all
available options:

.. click:: maud.cli:sample
   :prog: sample
   :show-nested:


Maud also has a command called :code:`simulate`, which is used in much the same
way as :code:`sample`:

.. click:: console

          maud simulate path/to/myinput


Here is a the full documentation for :code:`simulate`:

.. click:: maud.cli:sample
   :prog: simulate
   :show-nested:
