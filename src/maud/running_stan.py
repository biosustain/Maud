"""Code for sampling from a posterior distribution."""

import os

import arviz as az
import cmdstanpy
from cmdstanpy.stanfit.mcmc import CmdStanMCMC
from cmdstanpy.stanfit.vb import CmdStanVB

from maud.data_model.maud_input import MaudInput


HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_PRIOR_LOC_DRAIN = None
DEFAULT_PRIOR_SCALE_DRAIN = None
STAN_PROGRAM_RELATIVE_PATH = os.path.join("stan", "model.stan")
STAN_PROGRAM_RELATIVE_PATH_PREDICT = os.path.join(
    "stan", "out_of_sample_model.stan"
)
PPC_PROGRAM_RELATIVE_PATH = "out_of_sample_model.stan"

DEFAULT_SAMPLE_CONFIG = {
    "iter_warmup": 5,
    "iter_sampling": 5,
    "chains": 2,
    "max_treedepth": 11,
    "show_progress": True,
    "step_size": 0.025,
    "adapt_delta": 0.99,
    "save_warmup": True,
    "threads_per_chain": 1,
}
DEFAULT_VARIATIONAL_CONFIG = {
    "algorithm": "meanfield",
    "output_samples": 10,
    "require_converged": True,
}
SIM_CONFIG = {
    "chains": 1,
    "fixed_param": True,
    "iter_warmup": 0,
    "show_progress": False,
    "threads_per_chain": 1,
}


def sample(mi: MaudInput, output_dir: str) -> CmdStanMCMC:
    """Sample from the posterior defined by mi.

    :param mi: a MaudInput object
    :param output_dir: a string specifying where to save the output.
    """
    model = cmdstanpy.CmdStanModel(
        stan_file=os.path.join(HERE, STAN_PROGRAM_RELATIVE_PATH),
        cpp_options=mi.config.cpp_options,
        stanc_options=mi.config.stanc_options,
    )
    set_up_output_dir(output_dir, mi)
    sample_args: dict = {
        "data": os.path.join(output_dir, "input_data_train.json"),
        "inits": os.path.join(output_dir, "inits.json"),
        "output_dir": output_dir,
    }
    sample_args = {**sample_args, **DEFAULT_SAMPLE_CONFIG}
    if mi.config.cmdstanpy_config is not None:
        sample_args = {**sample_args, **mi.config.cmdstanpy_config}
    return model.sample(**sample_args)


def variational(mi: MaudInput, output_dir: str) -> CmdStanVB:
    """Do variational inference for the posterior defined by mi.

    :param mi: a MaudInput object
    :param output_dir: a string specifying where to save the output.
    """
    mi_options = (
        {}
        if mi.config.variational_options is None
        else mi.config.variational_options
    )
    model = cmdstanpy.CmdStanModel(
        stan_file=os.path.join(HERE, STAN_PROGRAM_RELATIVE_PATH),
        cpp_options=mi.config.cpp_options,
        stanc_options=mi.config.stanc_options,
    )
    set_up_output_dir(output_dir, mi)
    return model.variational(
        data=os.path.join(output_dir, "input_data_train.json"),
        inits=os.path.join(output_dir, "inits.json"),
        **{
            **DEFAULT_VARIATIONAL_CONFIG,
            **mi_options,
            **{"output_dir": output_dir},
        },
    )


def simulate(mi: MaudInput, output_dir: str, n: int) -> CmdStanMCMC:
    """Generate simulations from the prior mean.

    :param mi: a MaudInput object
    :param output_dir: a string specifying where to save the output.
    """
    model = cmdstanpy.CmdStanModel(
        stan_file=os.path.join(HERE, STAN_PROGRAM_RELATIVE_PATH),
        cpp_options=mi.config.cpp_options,
        stanc_options=mi.config.stanc_options,
    )
    set_up_output_dir(output_dir, mi)
    return model.sample(
        output_dir=output_dir,
        iter_sampling=n,
        data=os.path.join(output_dir, "input_data_train.json"),
        inits=os.path.join(output_dir, "inits.json"),
        **SIM_CONFIG,
    )


def set_up_output_dir(output_dir: str, mi: MaudInput):
    """Write input data and inits to the output directory."""
    input_path_train = os.path.join(output_dir, "input_data_train.json")
    input_path_test = os.path.join(output_dir, "input_data_test.json")
    inits_path = os.path.join(output_dir, "inits.json")
    cmdstanpy.utils.write_stan_json(input_path_train, mi.stan_input_train)
    if mi.stan_input_test is not None:
        cmdstanpy.utils.write_stan_json(input_path_test, mi.stan_input_test)
    cmdstanpy.utils.write_stan_json(inits_path, mi.inits_dict)


def predict(
    mi: MaudInput,
    output_dir: str,
    idata_train: az.InferenceData,
) -> az.InferenceData:
    """Call CmdStanModel.sample for out of sample predictions.

    :param mi: a MaudInput object
    :param output_dir: directory where output will be saved
    :param idata_train: InferenceData object with posterior draws
    """
    model = cmdstanpy.CmdStanModel(
        stan_file=os.path.join(HERE, STAN_PROGRAM_RELATIVE_PATH_PREDICT),
        cpp_options=mi.config.cpp_options,
        stanc_options=mi.config.stanc_options,
    )
    set_up_output_dir(output_dir, mi)
    kinetic_parameters = [
        "keq",
        "km",
        "kcat",
        "dissociation_constant",
        "transfer_constant",
        "kcat_phos",
        "ki",
    ]
    posterior = idata_train.get("posterior")
    sample_stats = idata_train.get("sample_stats")
    assert posterior is not None
    assert sample_stats is not None
    chains = sample_stats["chain"]
    draws = sample_stats["draw"]
    dims = {
        "conc_test": ["experiment", "mic"],
        "conc_enzyme_test": ["experiment", "enzyme"],
        "flux_test": ["experiment", "reaction"],
    }
    for chain in chains:
        for draw in draws:
            inits = {
                par: (
                    posterior[par]
                    .sel(chain=chain, draw=draw)
                    .to_series()
                    .values
                )
                for par in kinetic_parameters
                if par in posterior.keys()
            }
            sample_args: dict = {
                "data": os.path.join(output_dir, "input_data_test.json"),
                "inits": inits,
                "output_dir": output_dir,
                "iter_warmup": 0,
                "iter_sampling": 1,
                "fixed_param": True,
                "show_progress": False,
            }
            if mi.config.cmdstanpy_config_predict is not None:
                sample_args = {
                    **sample_args,
                    **mi.config.cmdstanpy_config_predict,
                }
            mcmc_draw = model.sample(**sample_args)
            idata_draw = az.from_cmdstan(
                mcmc_draw.runset.csv_files,
                coords={
                    "experiment": [e.id for e in mi.experiments if e.is_test],
                    "mic": [m.id for m in mi.kinetic_model.mics],
                    "enzyme": [e.id for e in mi.kinetic_model.enzymes],
                    "reaction": [r.id for r in mi.kinetic_model.reactions],
                },
                dims=dims,
            ).assign_coords(
                coords={"chain": [chain], "draw": [draw]},
                groups="posterior_groups",
            )
            if draw == 0:
                idata_chain = idata_draw.copy()
            else:
                idata_chain = az.concat(
                    [idata_chain, idata_draw], dim="draw", reset_dim=False
                )
        if chain == 0:
            out = idata_chain.copy()
        else:
            out = az.concat([out, idata_chain], dim="chain", reset_dim=False)
    return out
