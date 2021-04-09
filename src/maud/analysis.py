"""Functions for analysing Maud output."""

from math import ceil
from typing import Dict, List, Union

import arviz as az
import numpy as np
from matplotlib import pyplot as plt

from maud.data_model import MaudInput


def load_infd(csvs: List[str], mi: MaudInput) -> az.InferenceData:
    """Get an arviz InferenceData object from Maud csvs."""

    def join_list_of_strings(l1, l2, sep="-"):
        return list(map(lambda a: f"{a[0]}{sep}{a[1]}", zip(l1, l2)))

    coords = {
        **mi.stan_coords.__dict__,
        **{
            "kms": join_list_of_strings(mi.stan_coords.km_enzs, mi.stan_coords.km_mics),
            "yconcs": join_list_of_strings(
                mi.stan_coords.yconc_exps, mi.stan_coords.yconc_mics
            ),
            "yfluxs": join_list_of_strings(
                mi.stan_coords.yflux_exps, mi.stan_coords.yflux_rxns
            ),
            "yenzs": join_list_of_strings(
                mi.stan_coords.yenz_exps, mi.stan_coords.yenz_enzs
            ),
        },
    }
    return az.from_cmdstan(
        csvs,
        coords=coords,
        dims={
            "enzyme": ["experiments", "enzymes"],
            "conc": ["experiments", "mics"],
            "flux": ["experiments", "reactions"],
            "formation_energy": ["metabolites"],
            "kcat": ["enzymes"],
            "km": ["kms"],
            "yconc_sim": ["yconcs"],
            "yflux_sim": ["yfluxs"],
            "yenz_sim": ["yenzs"],
            "log_lik_conc": ["yconcs"],
            "log_lik_flux": ["yfluxs"],
            "log_lik_enz": ["yenzs"],
        },
        save_warmup=True,
    )


def plot_1d_var(
    infd: az.InferenceData,
    varname: str,
    codes: List[str],
    true_values: Union[List[float], None] = None,
    logscale: bool = False,
) -> tuple:
    """Plot a 1 dimensional variable."""
    tv = dict(zip(codes, true_values))
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
