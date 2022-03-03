"""Fetch multivariate formation energy priors from equilibrator.

NB: This script depends on the python packages equilibrator_api and
equilibrator_cache: make sure you have those installed in your python
environment!

"""

import argparse
import os
from typing import Tuple

import numpy as np
import pandas as pd
from equilibrator_api import ComponentContribution
from equilibrator_cache.models.compound import Compound

from maud.data_model import MaudInput
from maud.io import load_maud_input


HELP_MSG = """

    Use equilibrator to get a mean vector and covariance matrix for
    the formation energies of the metabolites of the kinetic model defined at
    maud_input_dir.

    Make sure that the metabolites-in-compartments in the kinetic model file
    have values in the field "metabolite_external_id" that equilibrator can
    recognise - otherwise this script will raise an error.

    NB: This script depends on the python packages equilibrator_api and
    equilibrator_cache: make sure you have those installed in your python
    environment!

"""


def get_dgf_priors(mi: MaudInput) -> Tuple[pd.Series, pd.DataFrame]:
    """Given a Maud input, get a multivariate prior from equilibrator.

    Returns a pandas Series of prior means and a pandas DataFrame of
    covariances. Both are indexed by metabolite ids.

    :param mi: A MaudInput object

    """
    cc = ComponentContribution()
    mu = []
    sigmas_fin = []
    sigmas_inf = []
    met_ix = pd.Index(mi.stan_coords.metabolites, name="metabolite")
    met_order = [m.id for m in mi.kinetic_model.metabolites]
    for m in mi.kinetic_model.metabolites:
        external_id = m.id if m.inchi_key is None else m.inchi_key
        c = cc.get_compound(external_id)
        if isinstance(c, Compound):
            mu_c, sigma_fin_c, sigma_inf_c = cc.standard_dg_formation(c)
            mu_c += c.transform(
                cc.p_h, cc.ionic_strength, cc.temperature, cc.p_mg
            ).m_as("kJ/mol")
            mu.append(mu_c)
            sigmas_fin.append(sigma_fin_c)
            sigmas_inf.append(sigma_inf_c)
        else:
            raise ValueError(
                f"cannot find compound for metabolite {m.id}"
                f" with external id {external_id}."
                "\nConsider setting the field metabolite_inchi_key"
                " if you haven't already."
            )
    sigmas_fin = np.array(sigmas_fin)
    sigmas_inf = np.array(sigmas_inf)
    cov = sigmas_fin @ sigmas_fin.T + 1e6 * sigmas_inf @ sigmas_inf.T
    cov = (
        pd.DataFrame(cov, index=met_order, columns=met_order)
        .loc[met_ix, met_ix]
        .round(10)
    )
    mu = pd.Series(mu, index=met_order, name="prior_mean_dgf").loc[met_ix].round(10)
    return mu, cov


def main():
    """Run the script."""
    parser = argparse.ArgumentParser(description=HELP_MSG)
    parser.add_argument(
        "maud_input_dir",
        type=str,
        nargs=1,
        help="A path to a Maud input directory",
    )
    maud_input_dir = parser.parse_args().maud_input_dir[0]
    if os.path.exists(maud_input_dir):
        mi = load_maud_input(maud_input_dir, "sample")
        mu, cov = get_dgf_priors(mi)
        print("Prior mean vector:")
        print(mu)
        print("Prior covariance:")
        print(cov)
        mu.to_csv(os.path.join(maud_input_dir, "dgf_prior_mean_equilibrator.csv"))
        cov.to_csv(os.path.join(maud_input_dir, "dgf_prior_cov_equilibrator.csv"))
    else:
        raise ValueError(f"{maud_input_dir} is not a valid path.")


if __name__ == "__main__":
    main()
