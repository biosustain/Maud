# 15. Change input format to a folder

date: 2021-01-20

## status

accepted

## context

Maud's input is currently a single [toml](https://github.com/toml-lang/toml)
file. 

In the long run we hope to generate input data for Maud automatically from 3D
genome scale models, raw measurement results and preliminary analyses to
generate priors, with minimal manual input to trim the model down to a
manageable size.

The toml format is fine for us to work with at the moment but has a few
drawbacks:

- Expressing prior information and experiment results in toml results in an
  unnecessarily complicated tree structure and also very long files.
- Our toml specification amounts to a custom language for specifying kinetic
  models and the required data: this is hard for new users to learn and limits
  Maud's interoperability.
- Toml files probably won't be used in the long term end state.

Recently Nick found a similar format called
[PEtab](https://github.com/PEtab-dev/PEtab) which seems to solve some of these
problems.


## decision

Move to a new format that is like PEtab in that kinetic models, experiment
results and priors are contained in different files and formats, but without
going all the way to inputting valid PEtab folders.

Specifically, under the proposal a Maud input directory must contain a file
called `config.toml` which contains paths to:

- a toml file representing a kinetic model with no change to the current
  format.
- a csv file with priors
- a csv file with measurements

## consequences

Good consequences:

- It will be easier to change the kinetic model format (or priors/experiments)
  without needing to change the other parts.
- All run-specific configuration can be expressed in files rather than as
  command line arguments. This is good for reproducibility and will also
  simplify the cli.
- We will be able to manipulate priors and experiments files with tools that
  are designed for csv files.

Bad consequences:

- We will probably need to change the input again when we integrate Maud with
  the CfB's data and standard format genome scale models.

- We will need to port models that are not currently in the Maud repository to
  the new format.

- We may need to update maudeller.
