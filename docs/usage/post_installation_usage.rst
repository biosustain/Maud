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

Maud also has a commands called :code:`simulate`, which is used in much the
same way as :code:`sample`. Here is the full documentation for both:

.. click:: maud.cli:sample_command
   :show-nested:
   :prog: sample


.. click:: maud.cli:simulate_command
   :show-nested:
   :prog: simulate
