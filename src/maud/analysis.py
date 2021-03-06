"""Functions for analysing Maud output."""

from math import ceil
from typing import Dict, List, Union

import arviz as az
import numpy as np
from matplotlib import pyplot as plt

from maud.data_model import MaudInput


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
        },
    )


def plot_1d_var(
    infd: az.InferenceData,
    varname: str,
    codes: dict,
    true_values: Union[List[float], None] = None,
    logscale: bool = False,
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
        x = samples[col]
        hist, bins = np.histogram(x, bins=30)
        xscale = "linear"
        if logscale:
            bins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
            xscale = "log"
        _, _, hist_patches = ax.hist(x, bins=bins)
        ax.set_xscale(xscale)
        vline = ax.axvline(tv[col], color="red")
        ax.set_title(col)
    f.legend(
        [hist_patches[0], vline],
        ["Marginal posterior", "True value"],
        frameon=False,
        loc="lower center",
        ncol=2,
    )
    return f, axes


def plot_experiment_var(
    infd: az.InferenceData,
    varname: str,
    var_codes: dict,
    exp_codes: dict,
    true_values: Union[Dict[str, Dict[str, float]], None] = None,
    meas_values: Union[Dict[str, Dict[str, float]], None] = None,
    logscale: bool = False,
) -> tuple:
    """Plot a 2d variable where the first dimension is experiment."""
    samples = infd.posterior[varname].to_series().unstack([-1, -2])
    var_labs, exp_labs = samples.columns.levels
    f, axes = plt.subplots(len(exp_labs), len(var_labs))
    for exp, axrow in zip(exp_labs, axes):
        for var, ax in zip(var_labs, axrow):
            x = samples[(var, exp)]
            hist, bins = np.histogram(x, bins=30)
            xscale = "linear"
            if logscale:
                bins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
                xscale = "log"
            _, _, hist_patches = ax.hist(x, bins=bins)
            ax.set_xscale(xscale)
            if true_values is not None:
                if exp in true_values.keys():
                    if var in true_values[exp].keys():
                        vline_truth = ax.axvline(true_values[exp][var], color="red")
            if meas_values is not None:
                if exp in meas_values.keys():
                    if var in meas_values[exp].keys():
                        vline_meas = ax.axvline(meas_values[exp][var], color="orange")
            ax.set_title(var)
        axrow[0].set_ylabel(exp)
    leg_handles = [hist_patches[0]]
    leg_labs = ["Marginal posterior"]
    if true_values is not None:
        leg_handles += [vline_truth]
        leg_labs += ["True value"]
    if meas_values is not None:
        leg_handles += [vline_meas]
        leg_labs += ["Measured value"]
    f.legend(leg_handles, leg_labs, frameon=False, loc="lower center", ncol=3)
    return f, axes
