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

==========================
Out-of-sample predictions
==========================

After generating posterior draws you can use these draws to test on future
experiments.


To generate predictions from testing data at :code:`path/to/my/testing/input`
using samples from training data at :code:`path/to/my/training/input` run the
following command:


.. code:: console

   maud generate-predictions path/to/my/training/input --oos_path="path/to/my/testing/input"


This will trigger maud to generate an output folder in the directory where the
command is run. To save the output folder somewhere else, use the option
:code:`output_dir` as follows:

.. code:: console

   <previous cmd> --output_dir path/to/desired/output/dir

Here is the full documentation:


.. click:: maud.cli:generate_predictions_command
   :show-nested:
   :prog: generate-predictions

