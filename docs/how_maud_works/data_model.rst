Maud's data model
=================

The main job of Maud's data model is to define an object called
:code:`MaudInput` which represents in Python all the information provided by the
user, and holds the information that is passed to Stan. It extensively uses
`pydantic dataclasses <https://pydantic-docs.helpmanual.io/usage/dataclasses/>`_
for validation.

Code implementing Maud's data model is in the folder 
`src/maud/data_model <https://github.com/biosustain/Maud/tree/master/src/maud/data_model>`_.

The information provided by the users consists of:

- Configuration represented as a :code:`MaudConfig` object
- A kinetic model represented by a :code:`KineticModel` object
- Priors represented by a :code:`UserPriorInput` object
- Measurements represented by a :code:`MeasurementSet` object
- Initial parameter values represented by an optional pandas :code:`DataFrame`

On initialisation a :code:`MaudInput` performs the following steps:

- Stan variables are inferred from the kinetic model and measurements, creating
  a :code:`StanVariableSet`.
- User-specified priors are supplemented with defaults to create a :code:`PriorSet`. The :code:`StanVariableSet` is required at this stage in order to ensure correct shapes.
- A full set of initial values is created using prior means as defaults.
- Stan inputs for the train and test models are created as ``StanInputTrain`` and ``StanInputTest`` objects.

`MaudConfig
<https://github.com/biosustain/Maud/tree/master/src/maud/data_model/maud_config.py>`_
is a simple python representation for the user toml input and specifies defaults
for the non-obligatory fields.

The `KineticModel
<https://github.com/biosustain/Maud/tree/master/src/maud/data_model/kinetic_model.py>`_
class validates the user toml input and supplements it with attributes
:code:`drains`, :code:`edges` and :code:`stoichiometric_matrix`.

`UserPriorInput
<https://github.com/biosustain/Maud/tree/master/src/maud/data_model/prior_set.py>`_
is a thin container for the user-provided prior files.

`MeasurementSet
<https://github.com/biosustain/Maud/tree/master/src/maud/data_model/prior_set.py>`_
contains user-provided measurement and experimental setup information.

`StanVariableSet
<https://github.com/biosustain/Maud/tree/master/src/maud/data_model/stan_variable_set.py>`_
specifies which parameters appear in Maud's statistical model, and provides
structural information about each parameter, such as its shape, default location
and scale and composition of its ids.

`PriorSet
<https://github.com/biosustain/Maud/tree/master/src/maud/data_model/prior_set.py>`_
specifies what kind of prior each parameter should have. The options are
independent 1d, independent 2d or multivariate 1d.

`StanInputTrain and StanInputTest
<https://github.com/biosustain/Maud/tree/master/src/maud/data_model/stan_input.py>`_
specify what data can possibly go into Maud's Stan input dictionaries, performs
validation and sets out what these dictionaries should look like.
