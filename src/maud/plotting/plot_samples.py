# Copyright (C) 2020 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Code for plotting posterior distribution."""
import argparse
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotnine as p9
import seaborn as sns

from maud import io
from maud.analysis import load_infd
from maud.data_model import IndPrior1d, IndPrior2d
from maud.user_templates import get_parameter_coords


HELP_MSG = """
This script plots violin plots of all parameters.
If there were priors for these parameters it will
plot them with the 95% CI, except for multivariate
parameters. If there are measurements they will be
plotted as a single point. This script also generates
pair plots of parameters within an enzyme.
"""

VARIABLES_TO_ANALYSE = [
    "kcat",
    "kcat_phos",
    "km",
    "drain",
    "ki",
    "diss_t",
    "diss_r",
    "transfer_constant",
    "conc_unbalanced",
    "conc_enzyme",
    "conc_phos",
    "conc",
    "flux",
    "keq",
    "dgf",
]

ENZYME_GROUP = [
    "kcat",
    "keq",
    "km",
    "ki",
    "transfer_constant",
    "diss_t",
    "diss_r",
]

LOG_SCALE_VARIABLES = [
    "kcat",
    "kcat_phos",
    "km",
    "ki",
    "keq",
    "diss_t",
    "diss_r",
    "transfer_constant",
    "conc_unbalanced",
    "conc_enzyme",
    "conc_phos",
    "conc",
]
UNITS = {
    "kcat": "1/s",
    "kcat_phos": "1/s",
    "km": "mM",
    "drain": "mM/s",
    "ki": "mM",
    "diss_t": "mM",
    "diss_r": "mM",
    "transfer_constant": "",
    "conc_unbalanced": "mM",
    "conc_enzyme": "mM",
    "conc_phos": "mM",
    "conc": "mM",
    "flux": "mM/s",
    "keq": "",
    "dgf": "kJ/mmol",
}


def get_dims_enz(par, parameter_coords, var_to_dims):
    """Return dataframe with enzyme_id and parameter_ids."""
    par_dataframe = pd.DataFrame()
    for p in parameter_coords:
        if p.id == par:
            if p.linking_list:
                par_dataframe["par_id"] = list(p.linking_list.values())[0]
                par_dataframe["enzyme_id"] = list(p.coords["enzyme_id"])
            else:
                if "enzyme_id" in p.coords.keys():
                    par_dataframe["par_id"] = list(p.coords["enzyme_id"])
                    par_dataframe["enzyme_id"] = list(p.coords["enzyme_id"])
                elif "edges" in p.coords.keys():
                    par_dataframe["par_id"] = list(p.coords["edges"])
                    par_dataframe["enzyme_id"] = list(p.coords["edges"])
    return par_dataframe


def get_ci_1d(p):
    """Return lower and upper CI given 1d prior."""
    if p.parameter_name in LOG_SCALE_VARIABLES:
        mean = np.log(p.location.reset_index()["location"].values)
        std = p.scale.reset_index()["scale"].values
        lower_ci = np.exp(mean - 2 * std)
        upper_ci = np.exp(mean + 2 * std)
    else:
        mean = p.location.reset_index()["location"].values
        std = p.scale.reset_index()["scale"].values
        lower_ci = mean - 2 * std
        upper_ci = mean + 2 * std

    return lower_ci, upper_ci


def plot_violin_plots(
    par_id: str,
    dims: List[str],
    draws: Dict,
    log_scale_variables: List[str],
    units: Dict[str, str],
    confidence_intervals,
    measurements,
):
    """Plot and save violin plots of parsed distributions.

    :param par_id: Name of the parameter plotted
    :param dims: Dimensions of the parameter
    :param draws: pd.Dataframe of parameter distribution
    indexed by dimensions and contains the population samples
    :param log_scale_variables: Parameters that are log-distributed
    :param units: Dictionary of units for each parameter
    """
    par_units = units[par_id]
    x = fill = dims[0] if len(dims) <= 1 else "experiments"
    plot = (
        p9.ggplot(data=draws)
        + p9.geom_violin(
            p9.aes(y=f"{par_id}", x=x, fill=fill),
            position="identity",
            color="None",
            size=0.5,
            alpha=0.7,
            weight=0.7,
            linetype="None",
        )
        + p9.labels.ylab(f"{par_id} {par_units}")
    )
    if par_id in confidence_intervals.keys():
        plot += p9.geoms.geom_errorbar(
            p9.aes(x=x, ymin="lower_ci", ymax="upper_ci"),
            data=confidence_intervals[par_id],
            width=0.1,
        )
    if par_id in measurements.keys():
        if len(measurements[par_id]) > 0:
            plot += p9.geoms.geom_point(
                p9.aes(y="measurement", x=x),
                data=measurements[par_id],
            )
    if len(dims) == 1:
        plot += p9.themes.theme(
            axis_text_x=p9.element_text(angle=70),
        )
    if len(dims) > 1:
        plot += p9.facet_wrap(f"~{dims[1]}") + p9.themes.theme(
            panel_spacing_y=0.05,
            panel_spacing_x=0.35,
            axis_title=p9.element_text(size=10),
            axis_text=p9.element_text(size=11),
            axis_text_y=p9.element_text(size=8, angle=45),
            axis_title_x=p9.element_blank(),
            axis_text_x=p9.element_blank(),
        )
    if par_id in log_scale_variables:
        plot += p9.scale_y_log10()

    return plot


def plot_posteriors(maud_output_dir, output_dir):
    """Plot posterior distributions of Maud model."""
    # Collecting information from draws and maud input
    csvs = list(Path(maud_output_dir / "samples").rglob("*.csv"))
    mi = io.load_maud_input(data_path=maud_output_dir / "user_input", mode="sample")
    parameter_coords = get_parameter_coords(mi.stan_coords)
    infd = load_infd(csvs, mi)
    list_of_model_variables = list(infd.posterior.variables.keys())
    var_to_dims = {
        var: list(infd.posterior[var].dims[2:])
        for var in VARIABLES_TO_ANALYSE
        if var in list_of_model_variables
    }
    var_to_draws = {
        var: infd.posterior[var].to_dataframe().reset_index()
        for var in VARIABLES_TO_ANALYSE
        if var in list_of_model_variables
    }
    enzyme_dims = {
        par: get_dims_enz(par, parameter_coords, var_to_dims)
        for par in ENZYME_GROUP
        if par in list_of_model_variables
    }
    priors = mi.priors
    confidence_intervals = dict()
    measurements = dict()
    # Retriving priors with confidence intervals (CIs)
    for par in parameter_coords:
        if par.id in list_of_model_variables:
            if f"priors_{par.id}" in dir(priors):
                par_dataframe = pd.DataFrame.from_dict(par.coords)
                if par.linking_list is None:
                    coords_rename = {
                        scs: infd_coord
                        for scs, infd_coord in zip(
                            list(par.coords.keys()), par.infd_coord_list
                        )
                    }
                    par_dataframe = par_dataframe.rename(columns=(coords_rename))
                else:
                    par_dataframe[par.infd_coord_list[0]] = list(
                        par.linking_list.values()
                    )[0]
                par_dataframe["parameter_name"] = par.id
                p = getattr(priors, f"priors_{par.id}")
                if isinstance(p, IndPrior1d):
                    lower_ci, upper_ci = get_ci_1d(p)
                    par_dataframe["lower_ci"] = lower_ci
                    par_dataframe["upper_ci"] = upper_ci
                    confidence_intervals[par.id] = par_dataframe
                elif isinstance(p, IndPrior2d):
                    location_df = (
                        p.location.unstack()
                        .reset_index()
                        .rename(columns=({0: "location"}))
                    )
                    scale_df = (
                        p.scale.unstack().reset_index().rename(columns=({0: "scale"}))
                    )
                    par_dataframe = par_dataframe.merge(
                        location_df,
                        left_on=par.infd_coord_list,
                        right_on=list(par.coords.keys()),
                    )
                    par_dataframe = par_dataframe.drop(list(par.coords.keys()), axis=1)
                    par_dataframe = par_dataframe.merge(
                        scale_df,
                        left_on=par.infd_coord_list,
                        right_on=list(par.coords.keys()),
                    )
                    par_dataframe = par_dataframe.drop(list(par.coords.keys()), axis=1)
                    if par.id in LOG_SCALE_VARIABLES:
                        par_dataframe["lower_ci"] = par_dataframe.apply(
                            lambda x: np.exp(np.log(x["location"]) - 2 * x["scale"]),
                            axis=1,
                        )
                        par_dataframe["upper_ci"] = par_dataframe.apply(
                            lambda x: np.exp(np.log(x["location"]) + 2 * x["scale"]),
                            axis=1,
                        )
                    else:
                        par_dataframe["lower_ci"] = par_dataframe.apply(
                            lambda x: x["location"] - 2 * x["scale"], axis=1
                        )
                        par_dataframe["upper_ci"] = par_dataframe.apply(
                            lambda x: x["location"] + 2 * x["scale"], axis=1
                        )
                    confidence_intervals[par.id] = par_dataframe
    # Retriving mean of measurement
    for measurement_id, measurement_type in zip(
        ["yconc", "yflux", "yenz"], ["conc", "flux", "conc_enzyme"]
    ):
        rename_columns = {"conc": "mics", "flux": "reactions", "conc_enzyme": "enzymes"}
        tmp_measurements = getattr(mi.measurements, measurement_id).reset_index()
        tmp_measurements = tmp_measurements.rename(
            columns=({"experiment_id": "experiments", "target_id": measurement_type})
        )
        tmp_measurements = tmp_measurements.rename(columns=(rename_columns))
        measurements[measurement_type] = tmp_measurements
    # Plotting violin plots from parameter distributions
    for var in list(var_to_dims.keys()):
        dims = var_to_dims[var]
        draws = var_to_draws[var]
        plot = plot_violin_plots(
            var,
            dims,
            draws,
            LOG_SCALE_VARIABLES,
            UNITS,
            confidence_intervals,
            measurements,
        )
        plot.save(
            filename=output_dir / f"{var}_posterior.png",
            verbose=False,
            dpi=300,
        )
    # plotting pairplots of enzyme parameters
    for enz in mi.stan_coords.enzymes:
        enz_par_df = pd.DataFrame()
        for par, par_df in enzyme_dims.items():
            par_draws = var_to_draws[par]
            enz_dims = par_df[par_df["enzyme_id"] == enz]["par_id"].to_list()
            if len(enz_dims) > 0:
                for par_ind in enz_dims:
                    tmp_enz_par_df = pd.DataFrame()
                    tmp_enz_par_df = par_draws.loc[
                        par_draws[var_to_dims[par][0]] == par_ind
                    ].copy()
                    enz_par_df[par + "-" + par_ind] = np.log(
                        tmp_enz_par_df[par].to_list()
                    )
        sns.pairplot(enz_par_df)
        plt.savefig(output_dir / f"{enz}_pairplot.png")


def main():
    """Run the script."""
    parser = argparse.ArgumentParser(description=HELP_MSG)
    parser.add_argument(
        "maud_output_dir", type=str, nargs=1, help="A path to a Maud output directory"
    )
    parser.add_argument("--output_dir", default=".", help="Path to output directory")
    args = parser.parse_args()
    maud_output_dir = Path.cwd() / args.maud_output_dir[0]
    output_dir = Path.cwd() / args.output_dir
    plot_posteriors(maud_output_dir, output_dir)


if __name__ == "__main__":
    main()
