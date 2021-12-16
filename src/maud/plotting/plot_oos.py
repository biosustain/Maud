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

import numpy as np
import pandas as pd
import plotnine as p9
from plotnine import aes, geom_point, ggplot, labs

from maud import io


HELP_MSG = """
This script plots boxplots of enzyme concentrations,
metabolite concentrations, and fluxes of out-of-sample
predictions. Measurements are also included in these
plots.
"""


def plot_box_plots(var, draws, measurements, variable_id_map):
    """Return plotnine.geoms.geom_boxplot of given variable."""
    plot = p9.ggplot(data=draws[var]) + p9.geom_boxplot(
        p9.aes(x=variable_id_map[var], y=var, fill=variable_id_map[var]),
        outlier_shape="",
    )
    if measurements[var].empty is False:
        plot += p9.geoms.geom_point(
            p9.aes(y="measurement", x=variable_id_map[var]), data=measurements[var]
        )
    if var != "flux":
        plot += p9.scale_y_log10()
    plot += p9.facet_wrap("~experiments") + p9.themes.theme(
        panel_spacing_y=0.05,
        panel_spacing_x=0.35,
        axis_title=p9.element_text(size=10),
        axis_text=p9.themes.element_text(size=11),
    )
    if var == "flux":
        plot += p9.scale_y_continuous(
            breaks=np.arange(-0.001, 0.002, 0.00025), limits=[-0.001, 0.002]
        )
    plot += p9.theme(axis_text_x=p9.themes.element_text(rotation=90, size=6))
    return plot


def plot_oos(maud_oos_dir, output_dir):
    """Save boxplots given maud oos predictions and ourput_dir."""
    mi = io.load_maud_input(data_path=maud_oos_dir / "user_input", mode="predict")
    flux_df = pd.read_csv(maud_oos_dir / "oos_samples" / "flux.csv")
    conc_df = pd.read_csv(maud_oos_dir / "oos_samples" / "conc.csv")
    conc_enzyme_df = pd.read_csv(maud_oos_dir / "oos_samples" / "conc_enzyme.csv")
    flux = (
        flux_df.groupby(["experiments", "reactions"])
        .mean()
        .reset_index()[["experiments", "reactions", "flux"]]
        .rename({"flux": "mean"}, axis=1)
    )
    flux["var"] = (
        flux_df.groupby(["experiments", "reactions"]).var().reset_index()["flux"]
    )
    flux["explore"] = np.abs(flux["var"] / flux["mean"])
    flux_exploration = (
        flux.groupby(["experiments"]).mean().reset_index()[["experiments", "explore"]]
    )
    exploration_plot = (
        ggplot(flux_exploration)
        + aes(x="experiments", y="explore")
        + geom_point(size=2)
        + labs(
            y="var / mean",
            x="experiment",
            title="Average exploration for each experiment condition",
        )
    )
    exploration_plot += p9.themes.theme(
        axis_text_x=p9.element_text(angle=70),
    )
    exploration_plot.save(
        filename="Exploration_list.png",
        verbose=False,
        dpi=300,
    )

    draws = {"flux": flux_df, "conc": conc_df, "conc_enzyme": conc_enzyme_df}
    variable_id_map = {"conc": "mics", "flux": "reactions", "conc_enzyme": "enzymes"}
    measurements = {}
    for measurement_id, measurement_type in zip(
        ["yconc", "yflux", "yenz"], ["conc", "flux", "conc_enzyme"]
    ):

        tmp_measurements = getattr(mi.measurements, measurement_id).reset_index()
        tmp_measurements = tmp_measurements.rename(
            columns=({"experiment_id": "experiments", "target_id": measurement_type})
        )
        tmp_measurements = tmp_measurements.rename(columns=(variable_id_map))
        measurements[measurement_type] = tmp_measurements

    for var in ["flux", "conc", "conc_enzyme"]:
        plot = plot_box_plots(var, draws, measurements, variable_id_map)
        plot.save(
            filename=output_dir / f"{var}_oos.png",
            verbose=False,
            dpi=300,
        )
    return


def main():
    """Run the script."""
    parser = argparse.ArgumentParser(description=HELP_MSG)
    parser.add_argument(
        "maud_oos_dir", type=str, nargs=1, help="A path to a Maud output directory"
    )
    parser.add_argument("--output_dir", default=".", help="Path to output directory")
    args = parser.parse_args()
    maud_oos_dir = Path.cwd() / args.maud_oos_dir[0]
    output_dir = Path.cwd() / args.output_dir
    plot_oos(maud_oos_dir, output_dir)


if __name__ == "__main__":
    main()
