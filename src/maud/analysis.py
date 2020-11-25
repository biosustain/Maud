"""Functions for analysing Maud output."""

import arviz as az
from typing import List, Union
from maud.data_model import MaudInput
import pandas as pd
from matplotlib import pyplot as plt
from math import ceil


def load_infd(csvs: List[str], mi: MaudInput) -> az.InferenceData:
    """Get an arviz InferenceData object from Maud csvs."""
    return az.from_cmdstan(
        csvs,
        coords={
            "enzyme_name": list(mi.stan_codes.enzyme_codes.keys()),
            "mic_name": list(mi.stan_codes.mic_codes.keys()),
            "reaction": list(mi.stan_codes.reaction_codes.keys()),
            "metabolite": list(mi.stan_codes.metabolite_codes.keys()),
            "experiment": list(mi.stan_codes.experiment_codes.keys()),
            "km_id": [p.id[3:] for p in mi.priors.km_priors],
        },
        dims={
            "enzyme": ["experiment", "enzyme_name"],
            "conc": ["experiment", "mic_name"],
            "flux": ["experiment", "reaction"],
            "formation_energy": ["metabolite"],
            "kcat": ["enzyme_name"],
            "km": ["km_id"],
        }
    ) 



def plot_1d_var(
    infd: az.InferenceData,
    varname: str,
    codes: dict,
    true_values: Union[List[float], None] = None
) -> tuple:
    """Plot a 1 dimensional variable."""
    tv = dict(zip(codes.keys(), true_values))
    samples = infd.posterior[varname].to_series().unstack()
    nrow = int(len(samples.columns) ** 0.5 // 1)
    ncol = ceil(len(samples.columns) / nrow)
    f, axes = plt.subplots(nrow, ncol)
    axes = axes.ravel()
    for i, col in enumerate(samples.columns):
        ax = axes[i]
        _, _, hist_patches = ax.hist(samples[col], bins=30)
        vline = ax.axvline(tv[col], color="red")
        ax.set_title(col)
    f.legend(
        [hist_patches[0], vline],["Marginal posterior", "True value"],
        frameon=False, loc="lower center", ncol=2
    )
    return f, axes
    


def plot_experiment_var(
    infd: az.InferenceData,
    varname: str,
    var_codes: dict,
    exp_codes: dict,
    true_values: Union[List[float], None] = None
) -> tuple:
    """Plot a 2d variable where the first dimension is experiment."""
    samples = infd.posterior[varname].to_series().unstack([-1, -2])
    var_labs, exp_labs = samples.columns.levels
    tvdf = pd.DataFrame(
        true_values, columns=var_codes.keys(), index=exp_codes.keys()
    )
    f, axes = plt.subplots(len(exp_labs), len(var_labs))
    for exp, axrow in zip(exp_labs, axes):
        for var, ax in zip(var_labs, axrow):
            _, _, hist_patches = ax.hist(samples[(var, exp)], bins=30)
            if var in var_codes.keys() and true_values is not None:
                vline = ax.axvline(tvdf.loc[exp, var], color="red")
            ax.set_title(var)
        axrow[0].set_ylabel(exp)
    if true_values is None:
        leg_handles = [hist_patches[0]]
        leg_labs = ["Marginal posterior"]
    else:
        leg_handles = [hist_patches[0], vline]
        leg_labs = ["Marginal posterior", "True value"]
    f.legend(leg_handles, leg_labs, frameon=False, loc="lower center", ncol=2)
    return f, axes

