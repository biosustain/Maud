"""Functions for analysing Maud output."""

from math import ceil
from typing import Dict, List, Union

import arviz as az
import numpy as np
from arviz.utils import Numba
from matplotlib import pyplot as plt

from maud.data_model.maud_input import MaudInput


def join_list_of_strings(l1, l2, sep="-"):
    """Join strings for use in infd coordinates."""
    return list(map(lambda a: f"{a[0]}{sep}{a[1]}", zip(l1, l2)))


def load_infd_fit(fit, mi: MaudInput) -> az.InferenceData:
    """Get an arviz InferenceData object from out-of-sample tmp generated csvs."""

    coords = {
        **mi.stan_coords.__dict__,
        **{
            "reactions": mi.stan_coords.reactions,
            "kms": join_list_of_strings(
                mi.stan_coords.km_enzs, mi.stan_coords.km_mics
            ),
            "kis": join_list_of_strings(
                mi.stan_coords.ci_enzs, mi.stan_coords.ci_mics
            ),
            "diss_ts": join_list_of_strings(
                mi.stan_coords.ai_enzs, mi.stan_coords.ai_mics
            ),
            "diss_rs": join_list_of_strings(
                mi.stan_coords.aa_enzs, mi.stan_coords.aa_mics
            ),
        },
    }
    return az.from_cmdstan(
        fit,
        coords=coords,
        dims={
            "flux": ["experiments", "reactions"],
            "conc": ["experiments", "mics"],
            "conc_enzyme": ["experiments", "enzymes"],
            "conc_unbalanced": ["experiments", "unbalanced_mics"],
            "conc_phos": ["experiments", "phos_enzs"],
            "saturation": ["experiments", "edges"],
            "allostery": ["experiments", "edges"],
            "phosphorylation": ["experiments", "edges"],
            "reversibility": ["experiments", "edges"],
        },
        save_warmup=True,
    )


def load_infd(csvs: List[str], mi: MaudInput) -> az.InferenceData:
    """Get an arviz InferenceData object from Maud csvs."""

    Numba.disable_numba()
    coords = {
        **mi.stan_coords.__dict__,
        **{
            "reactions": mi.stan_coords.reactions,
            "kms": join_list_of_strings(
                mi.stan_coords.km_enzs, mi.stan_coords.km_mics
            ),
            "kis": join_list_of_strings(
                mi.stan_coords.ci_enzs, mi.stan_coords.ci_mics
            ),
            "diss_ts": join_list_of_strings(
                mi.stan_coords.ai_enzs, mi.stan_coords.ai_mics
            ),
            "diss_rs": join_list_of_strings(
                mi.stan_coords.aa_enzs, mi.stan_coords.aa_mics
            ),
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
            "flux": ["experiments", "reactions"],
            "conc": ["experiments", "mics"],
            "conc_enzyme": ["experiments", "enzymes"],
            "conc_unbalanced": ["experiments", "unbalanced_mics"],
            "conc_phos": ["experiments", "phos_enzs"],
            "drain": ["experiments", "drains"],
            "diss_t": ["diss_ts"],
            "diss_r": ["diss_rs"],
            "transfer_constant": ["allosteric_enzymes"],
            "dgf": ["metabolites"],
            "dgrs": ["experiments", "edges"],
            "keq": ["experiments", "edges"],
            "kcat": ["enzymes"],
            "kcat_phos": ["phos_enzs"],
            "km": ["kms"],
            "ki": ["kis"],
            "yconc_sim": ["yconcs"],
            "yflux_sim": ["yfluxs"],
            "yenz_sim": ["yenzs"],
            "log_lik_conc": ["yconcs"],
            "log_lik_flux": ["yfluxs"],
            "log_lik_enz": ["yenzs"],
            "saturation": ["experiments", "edges"],
            "allostery": ["experiments", "edges"],
            "phosphorylation": ["experiments", "edges"],
            "reversibility": ["experiments", "edges"],
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
                bins = np.logspace(
                    np.log10(bins[0]), np.log10(bins[-1]), len(bins)
                )
                xscale = "log"
            _, _, hist_patches = ax.hist(x, bins=bins)
            ax.set_xscale(xscale)
            if true_values is not None:
                if exp in true_values.keys():
                    if var in true_values[exp].keys():
                        vline_truth = ax.axvline(
                            true_values[exp][var], color="red"
                        )
            if meas_values is not None:
                if exp in meas_values.keys():
                    if var in meas_values[exp].keys():
                        vline_meas = ax.axvline(
                            meas_values[exp][var], color="orange"
                        )
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
