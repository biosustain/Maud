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
import os
from typing import Dict, List

import arviz as az
import plotnine as p9

from maud import io


HERE = os.path.dirname(os.path.abspath(__file__))
HOME = os.path.join(HERE, "..")
model_name = ""  # Toml name of file
run_number = ""  # Number specific to each run, found between "-"
PATHS = {
    "data": os.path.join(
        HOME, ""
    ),  # Where the input toml file is stored, relative to the Maud folder
    "results": os.path.join(
        HOME, ""
    ),  # The output directory, relative toe the Maud folder
}
CMDSTAN_FILE_TEMPLATE = os.path.join(
    PATHS["results"], "inference_model-{run_number}-{i}.csv"
)
N_CHAINS = 4
VARIABLES_TO_ANALYSE = [
    "conc",
    "flux",
    "keq",
    "kcat",
    "formation_energy",
    "km",
    "enzyme",
]
LOG_SCALE_VARIABLES = ["conc", "keq", "kcat", "km", "enzyme"]
UNITS = {
    "conc": "mM",
    "flux": "mM/s",
    "keq": "",
    "kcat": "1/s",
    "formation_energy": "kJ/mol",
    "km": "mM",
    "enzyme": "mM",
}


def plot_violin_plots(
    par_id: str,
    dims: List[str],
    draws: Dict,
    log_scale_variables: List[str],
    units: Dict[str, str],
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
            size=0.5,
            alpha=0.7,
            weight=0.7,
            linetype="None",
        )
        + p9.labels.ylab(f"{par_id} {par_units}")
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
        plot = plot + p9.scale_y_log10()
    return plot


def main():
    """Plot posterior distributions of Maud model."""
    files = [
        CMDSTAN_FILE_TEMPLATE.format(run_number=run_number, i=str(i))
        for i in range(1, N_CHAINS + 1)
    ]
    mi = io.load_maud_input_from_toml(os.path.join(PATHS["data"], f"{model_name}.toml"))
    infd = az.from_cmdstan(
        files,
        coords={
            "mics": list(mi.stan_codes["metabolite_in_compartment"].keys()),
            "mets": list(mi.stan_codes["metabolite"].keys()),
            "kms": [f"{p.enzyme_id}_{p.mic_id}" for p in mi.priors["kms"]],
            "enzymes": list(mi.stan_codes["enzyme"].keys()),
            "reactions": list(mi.stan_codes["reaction"].keys()),
            "experiments": list(mi.stan_codes["experiment"].keys()),
        },
        dims={
            "conc": ["experiments", "mics"],
            "flux": ["experiments", "reactions"],
            "keq": ["enzymes"],
            "kcat": ["enzymes"],
            "formation_energy": ["mets"],
            "km": ["kms"],
            "enzyme": ["experiments", "enzymes"],
        },
    )
    var_to_dims = {
        var: list(infd.posterior[var].dims[2:]) for var in VARIABLES_TO_ANALYSE
    }
    var_to_draws = {
        var: infd.posterior[var].to_dataframe().reset_index()
        for var in VARIABLES_TO_ANALYSE
    }
    for var in VARIABLES_TO_ANALYSE:
        dims = var_to_dims[var]
        draws = var_to_draws[var]
        plot = plot_violin_plots(var, dims, draws, LOG_SCALE_VARIABLES, UNITS)
        plot.save(
            filename=os.path.join(PATHS["results"], f"{var}_posterior.png"),
            verbose=False,
            dpi=300,
        )


if __name__ == "__main__":
    main()
