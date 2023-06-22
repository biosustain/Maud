============
Running Maud
============

Maud provides the following commands:

- :code:`maud sample [PATH TO MAUD INPUT]`
- :code:`maud simulate [PATH TO MAUD INPUT]`
- :code:`maud variational [PATH TO MAUD INPUT]`
- :code:`maud optimize [PATH TO MAUD INPUT]`
- :code:`maud predict [PATH TO MAUD SAMPLE OUTPUT]`

The documentation for each of these commands is as follows:

.. click:: maud.cli:sample_command
   :show-nested:
   :prog: sample

.. click:: maud.cli:simulate_command
   :show-nested:
   :prog: simulate

.. click:: maud.cli:variational_command
   :show-nested:
   :prog: variational

.. click:: maud.cli:optimize_command
   :show-nested:
   :prog: optimize

.. click:: maud.cli:predict_command
   :show-nested:
   :prog: predict
