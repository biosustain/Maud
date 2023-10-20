"""Provides model MaudConfig."""

from typing import Optional

from pydantic import BaseModel, ConfigDict, Field


class ODESolverConfig(BaseModel):
    """Config that is specific to an ODE solver."""

    rel_tol: float = 1e-9
    abs_tol: float = 1e-9
    max_num_steps: int = int(1e7)
    model_config: ConfigDict = {"frozen": True}


class AlgebraSolverConfig(BaseModel):
    """Config that is specific to an ODE solver."""

    rel_tol: float = 1e-7
    abs_tol: float = 1e-7
    max_num_steps: int = int(1e6)
    model_config: ConfigDict = {"frozen": True}


class MaudConfig(BaseModel):
    """User's configuration for a Maud input.

    :param name: name for the input. Used to name the output directory
    :param kinetic_model_file: path to a valid kientic model file.
    :param priors_file: path to a valid priors file.
    :param experiments_file: path to a valid experiments file.
    :param likelihood: Whether or not to take measurements into account.
    :param cmdstanpy_config: Arguments to cmdstanpy.CmdStanModel.sample.
    :param penalize_non_steady: Penalize the deviation from steady state in the log likelihood.
    :param ode_solver_config: Configuration for Stan's ode solver.
    :param algebra_solver_config: Configuration for Stan's algebra solver.
    :param stanc_options: Options for CmdStanModel argument `stanc_options`.
    :param cpp_options: Options for CmdStanModel `cpp_options`.
    :param variational_options: Arguments for CmdStanModel.variational.
    :param pathfinder_options: Arguments for CmdStanModel.pathfinder.
    :param laplace_options: Arguments for CmdStanModel.laplace.
    :param optimize_options: Arguments for CmdStanModel.optimize.
    :param user_inits_file: path to a csv file of initial values.
    :param steady_state_threshold_abs: abs threshold for Sv=0 be at steady state
    :param steady_state_threshold_rel: rel threshold for Sv=0 be at steady state
    :param steady_state_penalty_rel: relative standard deviation for SV checks when penalize_non_steady is true.
    :param default_initial_concentration: in molecule_unit per volume_unit
    :param drain_small_conc_corrector: number for correcting small conc drains
    :param molecule_unit: A unit for counting molecules, like 'mol' or 'mmol'
    :param volume_unit: A unit for measuring volume, like 'L'
    :param energy_unit: A unit for measuring energy, like 'J' or 'kJ'
    """  # noqa: E501

    name: str
    kinetic_model_file: str
    priors_file: str
    experiments_file: str
    likelihood: bool
    cmdstanpy_config: Optional[dict] = None
    cmdstanpy_config_predict: Optional[dict] = None
    stanc_options: Optional[dict] = None
    cpp_options: Optional[dict] = None
    variational_options: Optional[dict] = None
    pathfinder_options: Optional[dict] = None
    laplace_options: Optional[dict] = None
    optimize_options: Optional[dict] = None
    user_inits_file: Optional[str] = None
    ode_solver_config: ODESolverConfig = Field(default_factory=ODESolverConfig)
    algebra_solver_config: AlgebraSolverConfig = Field(
        default_factory=AlgebraSolverConfig
    )
    penalize_non_steady: bool = False
    steady_state_threshold_abs: float = 1e-8
    steady_state_threshold_rel: float = 1e-3
    steady_state_penalty_rel: float = 1e-8
    default_initial_concentration: float = 0.01
    drain_small_conc_corrector: float = 1e-6
    molecule_unit: str = "mmol"
    volume_unit: str = "L"
