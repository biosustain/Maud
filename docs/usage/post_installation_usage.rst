==========================
Generating posterior draws
==========================

Maud's main purpose is to generate an ensemble of samples that collectively
represent the information contained in some input data. 

To generate samples for input data at :code:`path/to/myinput.toml`, run the
following command:


.. code:: console
          maud sample path/to/data --n_chains=3 --n_warmup=200 --n_samples=100


This will use the data file at `path/to/myinput.toml` to create a Stan program
at `maud/stan_code/autogen/inference_model_myinput.stan`, compile it into a C++
Stan model, then start three independent MCMC chains, each using 200 iterations
for adaptation and storing the contents of 100 further iterations in a
file called `data/out/model_output_myinput-<n>.csv`.

Here is the full documentation for the :code:`sample` command, listing all
available options:

.. click:: maud.cli:sample
   :prog: sample
   :show-nested:

