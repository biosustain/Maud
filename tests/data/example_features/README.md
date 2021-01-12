# Feature Tests
This folder contains tests for the current mechanistic features
present in Maud. Using these defined models, `tests/test_linear.py`
runs an integration tests asserting that changes to maud stil result
in the recovery of the "true parameters." The true parameters can vary
as the measurement values are independent of the `toml` input file, with
the dependent parameters (`balanced_conc` and `flux`) being calculated
as a result of the input file values.

# List of features
Below is a list of features currently checked:

* Michaelis-Menten
* Allostery
	* Activating
	* Inhibiting
* Drains